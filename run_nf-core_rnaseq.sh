#!/usr/bin/env bash
#SBATCH --job-name=gtex_rnaseq_batches
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=logs/gtex_rnaseq_batches.%j.out
#SBATCH --error=logs/gtex_rnaseq_batches.%j.err

set -euo pipefail

SAMTOOLS_CPUS="${SLURM_CPUS_PER_TASK:-8}"
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

# Process each batch in this single SLURM allocation
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
    
    echo "Running nf-core/rnaseq for ${batch_id}"

    RIBO_INDEX=$(grep -oP '"ribo_database_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)
    RIBO_MANIFEST=$(grep -oP '"ribo_database_manifest"\s*:\s*"\K[^"]+' "conf/params.json" || true)
    CONTAM_INDEX=$(grep -oP '"contamination_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)

    NEXTFLOW_CONFIG_ARGS=()
    [ -f "conf/rnaseq.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/rnaseq.config" )
    [ -f "conf/star.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/star.config" )

    nextflow run rnaseq/main.nf \
        -params-file "conf/params.json" \
        "${NEXTFLOW_CONFIG_ARGS[@]}" \
        -profile singularity \
        --input "${samplesheet}" \
        --outdir "results/${batch_id}" \
        --contaminant_screening kraken2_bracken \
        --kraken_db "${CONTAM_INDEX}" \
        --without-wave \
        --save_unaligned \
        --skip_bbsplit \
        --skip_fastqc \
        --skip_rseqc \
        --skip_qualimap \
        --skip_dupradar \
        --skip_preseq \
        --skip_biotype_qc \
        --skip_stringtie \
        --skip_deseq2_qc \
        --skip_markduplicates \
        --skip_bigwig \
        --remove_ribo_rna \
        --ribo_removal_tool sortmerna \
        --ribo_database_manifest "${RIBO_MANIFEST}" \
        --sortmerna_index "${RIBO_INDEX}" \
        -with-trace "${batch_id}_trace.tsv" \
        -work-dir "work_nfcore/${batch_id}" \
        -with-tower \
        -name "nf_core_gtex_rnaseq_${batch_id}"

    echo "nf-core/rnaseq complete for ${batch_id}"
    echo "Cleaning up FASTQ files..."
    rm -rf "${batch_fastq_dir}"
    echo "Cleanup complete for ${batch_id}"

done

echo "All batches finished."
