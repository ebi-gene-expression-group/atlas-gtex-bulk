#!/usr/bin/env bash
set -euo pipefail

NFCORE_VERSION="3.18.0"
PROFILE="singularity,slurm"
GENOME="GRCh38"
ALIGNER="star_salmon"
SAMTOOLS_CPUS=8
SAMTOOLS_MEM="32G"
FASTQ_OUTDIR="fastq_batches"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-}"

if [ -z "${SAMTOOLS_IMAGE}" ]; then
    echo "ERROR: SAMTOOLS_IMAGE environment variable not set"
    exit 1
fi

mkdir -p logs results work_nfcore "${FASTQ_OUTDIR}"

# Function to convert BAMs to FASTQ for a batch
convert_batch_to_fastq() {
    local batch_info="$1"
    local batch_id=$(basename "${batch_info}" _info.csv)
    local batch_fastq_dir="${FASTQ_OUTDIR}/${batch_id}"
    
    mkdir -p "${batch_fastq_dir}"
    
    echo "Converting BAMs to FASTQ for ${batch_id}..."
    
    # Skip header row and process each BAM
    tail -n +2 "${batch_info}" | while IFS=',' read -r sample bam strandedness read_type; do
        echo "  Converting ${sample} (${read_type})..."

        singularity exec ${SAMTOOLS_IMAGE} samtools collate \
            -@ ${SAMTOOLS_CPUS} \
            -u \
            -O "${bam}" \
        | singularity exec ${SAMTOOLS_IMAGE} samtools fastq \
            -@ ${SAMTOOLS_CPUS} \
            -F 0x900 \
            -1 "${batch_fastq_dir}/${sample}_R1.fastq.gz" \
            -2 "${batch_fastq_dir}/${sample}_R2.fastq.gz" \
            -0 /dev/null \
            -s /dev/null \
            -n -
    done
    
    echo "FASTQ conversion complete for ${batch_id}"
}

# Function to create nf-core/rnaseq samplesheet from batch info and FASTQs
create_samplesheet() {
    local batch_info="$1"
    local batch_id=$(basename "${batch_info}" _info.csv)
    local batch_fastq_dir="${FASTQ_OUTDIR}/${batch_id}"
    local samplesheet="${batch_fastq_dir}/samplesheet.csv"
    
    echo "Creating samplesheet for ${batch_id}..."
    
    # Write header
    echo "sample,fastq_1,fastq_2,strandedness" > "${samplesheet}"
    
    # Process each BAM entry
    tail -n +2 "${batch_info}" | while IFS=',' read -r sample bam strandedness read_type; do
        echo "${sample},${batch_fastq_dir}/${sample}_R1.fastq.gz,${batch_fastq_dir}/${sample}_R2.fastq.gz,${strandedness}" >> "${samplesheet}"
    done
    
    echo "Samplesheet: ${samplesheet}"
    head -n 5 "${samplesheet}"
}

# Process each batch
for batch_info in batch_info/batch_*.csv; do
    [ -f "${batch_info}" ] || continue
    
    batch_id=$(basename "${batch_info}" _info.csv)
    batch_fastq_dir="${FASTQ_OUTDIR}/${batch_id}"
    
    echo "=========================================="
    echo "Processing ${batch_id}"
    echo "=========================================="
    
    # Convert BAMs to FASTQ
    convert_batch_to_fastq "${batch_info}"
    
    # Create samplesheet
    samplesheet="${batch_fastq_dir}/samplesheet.csv"
    create_samplesheet "${batch_info}"
    
    # Submit nf-core/rnaseq job
    echo "Submitting nf-core/rnaseq for ${batch_id}"
    
    sbatch \
      --job-name="rnaseq_${batch_id}" \
      --cpus-per-task=1 \
      --mem=4G \
      --time=10-00:00:00 \
      --output="logs/${batch_id}.%j.out" \
      --error="logs/${batch_id}.%j.err" \
      <<EOF
#!/usr/bin/env bash
set -euo pipefail

BATCH_ID="${batch_id}"
BATCH_FASTQ_DIR="${batch_fastq_dir}"
SAMPLESHEET="${samplesheet}"

echo "Running nf-core/rnaseq for \${BATCH_ID}..."

nextflow run nf-core/rnaseq \\
  -r ${NFCORE_VERSION} \\
  -profile ${PROFILE} \\
  --input "\${SAMPLESHEET}" \\
  --outdir "results/\${BATCH_ID}" \\
  --genome ${GENOME} \\
  --aligner ${ALIGNER} \\
  -work-dir "work_nfcore/\${BATCH_ID}" \\
  -resume

echo "nf-core/rnaseq complete for \${BATCH_ID}"
echo "Cleaning up FASTQ files..."
rm -rf "\${BATCH_FASTQ_DIR}"
echo "Cleanup complete for \${BATCH_ID}"

EOF

done

echo "All batches submitted!"
