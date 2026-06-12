#!/usr/bin/env bash
#SBATCH --job-name=gtex_rnaseq_batches
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=logs/gtex_rnaseq_batches.%j.out
#SBATCH --error=logs/gtex_rnaseq_batches.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=anilthanki@ebi.ac.uk

set -euo pipefail

SAMTOOLS_CPUS="${SLURM_CPUS_PER_TASK:-8}"
FASTQ_OUTDIR="fastq_batches"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-}"

if [ -z "${SAMTOOLS_IMAGE}" ]; then
    echo "ERROR: SAMTOOLS_IMAGE environment variable not set"
    exit 1
fi

mkdir -p logs results work_nfcore "${FASTQ_OUTDIR}"

# Convert a single BAM to paired-end FASTQ
convert_one_bam() {
    local sample="$1"
    local bam="$2"
    local batch_fastq_dir="$3"
    local threads="$4"
    
    echo "    [$(date '+%Y-%m-%d %H:%M:%S')] Converting ${sample}..."
    singularity exec ${SAMTOOLS_IMAGE} samtools collate \
        -@ ${threads} \
        -u \
        -O "${bam}" \
    | singularity exec ${SAMTOOLS_IMAGE} samtools fastq \
        -@ ${threads} \
        -F 0x900 \
        -1 "${batch_fastq_dir}/${sample}_R1.fastq.gz" \
        -2 "${batch_fastq_dir}/${sample}_R2.fastq.gz" \
        -0 /dev/null \
        -s /dev/null \
        -n -
}

# Function to convert BAMs to FASTQ for a batch in parallel
convert_batch_to_fastq() {
    local batch_info="$1"
    local batch_id=$(basename "${batch_info}" _info.csv)
    local batch_fastq_dir="${FASTQ_OUTDIR}/${batch_id}"
    local MAX_PARALLEL=$(( (SAMTOOLS_CPUS + 1) / 2 ))  # ~2 threads per parallel job
    local threads=$(( SAMTOOLS_CPUS / MAX_PARALLEL ))
    
    mkdir -p "${batch_fastq_dir}"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converting BAMs to FASTQ for ${batch_id} (${MAX_PARALLEL} parallel jobs, ${threads} threads each)..."
    
    local pids=()
    local failed=0
    
    # Skip header row and spawn parallel conversion jobs
    tail -n +2 "${batch_info}" | while IFS=',' read -r sample bam strandedness read_type; do
        # Wait if we've reached max parallel jobs
        while [ $(jobs -r | wc -l) -ge ${MAX_PARALLEL} ]; do
            sleep 1
        done
        
        # Launch conversion in background
        convert_one_bam "${sample}" "${bam}" "${batch_fastq_dir}" "${threads}" &
    done
    
    # Wait for all background jobs to complete
    wait
    local wait_status=$?
    
    if [ ${wait_status} -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: FASTQ conversion failed for ${batch_id}"
        return 1
    fi
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FASTQ conversion complete for ${batch_id}"
    return 0
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

# Return success when all expected PE FASTQ files for a batch already exist.
batch_fastqs_exist() {
    local batch_info="$1"
    local batch_fastq_dir="$2"

    tail -n +2 "${batch_info}" | while IFS=',' read -r sample bam strandedness read_type; do
        if [ ! -s "${batch_fastq_dir}/${sample}_R1.fastq.gz" ] || [ ! -s "${batch_fastq_dir}/${sample}_R2.fastq.gz" ]; then
            exit 1
        fi
    done
}

# Process each batch in this single SLURM allocation
for batch_info in batch_info/batch_*.csv; do
    [ -f "${batch_info}" ] || continue
    
    batch_id=$(basename "${batch_info}" _info.csv)
    batch_fastq_dir="${FASTQ_OUTDIR}/${batch_id}"
    batch_done_marker="results/${batch_id}/.done"
    
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${batch_id}"
    echo "=========================================="
    
    if [ -f "${batch_done_marker}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch ${batch_id} already completed (marker found); skipping."
        continue
    fi
    
    # Convert BAMs to FASTQ (skip if samplesheet already exists)
    samplesheet="${batch_fastq_dir}/samplesheet.csv"

    if [ -s "${samplesheet}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Samplesheet already exists for ${batch_id}; skipping FASTQ conversion."
    else
        if batch_fastqs_exist "${batch_info}" "${batch_fastq_dir}"; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] All FASTQs already exist for ${batch_id}; skipping FASTQ conversion."
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting FASTQ conversion for ${batch_id}..."
            convert_batch_to_fastq "${batch_info}" || {
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: FASTQ conversion failed for ${batch_id}; skipping this batch."
                continue
            }
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] FASTQ conversion completed for ${batch_id}."
        fi
        create_samplesheet "${batch_info}"
    fi
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting nf-core/rnaseq for ${batch_id}"

    RIBO_INDEX=$(grep -oP '"ribo_database_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)
    RIBO_MANIFEST=$(grep -oP '"ribo_database_manifest"\s*:\s*"\K[^"]+' "conf/params.json" || true)
    CONTAM_INDEX=$(grep -oP '"contamination_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)

    NEXTFLOW_CONFIG_ARGS=()
    [ -f "conf/rnaseq.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/rnaseq.config" )
    [ -f "conf/star.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/star.config" )

    module load nextflow/25.04.6

    nextflow run rnaseq/main.nf \
        -params-file "conf/params.json" \
        "${NEXTFLOW_CONFIG_ARGS[@]}" \
        -profile singularity \
        --input "${samplesheet}" \
        --outdir "results/${batch_id}" \
	--igenomes_ignore \
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

    if [ $? -eq 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core/rnaseq completed successfully for ${batch_id}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleaning up FASTQ files and work directory..."
        rm -rf "${batch_fastq_dir}" "work_nfcore/${batch_id}" ".nextflow/history" || true
        mkdir -p "results/${batch_id}"
        touch "${batch_done_marker}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup complete for ${batch_id}; marked batch as complete."
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: nf-core/rnaseq failed for ${batch_id}; preserving FASTQ and work files for debugging."
    fi

done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All batches finished."
