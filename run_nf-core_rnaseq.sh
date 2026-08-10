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

FASTQ_OUTDIR="fastq_batches"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:?SAMTOOLS_IMAGE environment variable not set}"

mkdir -p logs results work_nfcore "${FASTQ_OUTDIR}"

declare -a nfcore_job_ids=()
declare -a nfcore_batches=()

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
    local sample=""
    local bam=""
    local strandedness=""
    local read_type=""

    while IFS=',' read -r sample bam strandedness read_type; do
        [ -z "${sample}" ] && continue
        if [ ! -s "${batch_fastq_dir}/${sample}_R1.fastq.gz" ] || [ ! -s "${batch_fastq_dir}/${sample}_R2.fastq.gz" ]; then
            return 1
        fi
    done < <(tail -n +2 "${batch_info}")

    return 0
}

# Process each batch: submit conversions, then start nf-core in background
for batch_info in batch_info/batch_*.csv; do
    [ -f "${batch_info}" ] || continue
    
    batch_id=$(basename "${batch_info}" _info.csv)
    batch_fastq_dir="${FASTQ_OUTDIR}/${batch_id}"
    batch_done_marker="results/${batch_id}/.done"
    batch_running_marker="results/${batch_id}/.running"
    batch_failed_marker="results/${batch_id}/.failed"
    samplesheet="${batch_fastq_dir}/samplesheet.csv"
    
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${batch_id}"
    echo "=========================================="
    
    if [ -f "${batch_done_marker}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch ${batch_id} already completed; skipping."
        continue
    fi

    if [ -f "${batch_running_marker}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch ${batch_id} has a running marker; skipping rerun."
        continue
    fi

    if [ -f "${batch_failed_marker}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch ${batch_id} has a failed marker; skipping rerun."
        continue
    fi
    
    mkdir -p "${batch_fastq_dir}"
    
    # Check if FASTQs already exist
    if [ -s "${samplesheet}" ] && batch_fastqs_exist "${batch_info}" "${batch_fastq_dir}"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FASTQs already exist for ${batch_id}; skipping conversion."
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting conversion jobs for ${batch_id} samples..."
        
        declare -a job_ids=()
        
        # Submit one job per sample in this batch
        while IFS=',' read -r sample bam strandedness read_type; do
            [ -z "${sample}" ] && continue
            
            job_id=$(sbatch --parsable \
                --export=ALL,SAMPLE="${sample}",BAM="${bam}",BATCH_FASTQ_DIR="${batch_fastq_dir}",SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE}" \
                convert_sample_fastq.sbatch)
            
            job_ids+=("${job_id}")
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitted job ${job_id} for sample ${sample}"
        done < <(tail -n +2 "${batch_info}")
        
        # Wait for all sample conversion jobs to complete
        if [ "${#job_ids[@]}" -gt 0 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Waiting for ${#job_ids[@]} sample conversion jobs to finish..."
            
            job_list=$(IFS=,; echo "${job_ids[*]}")
            while squeue -j "${job_list}" -h -o "%i" 2>/dev/null | grep -q .; do
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Still waiting..."
                sleep 10
            done
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] All sample conversion jobs completed for ${batch_id}"
        fi
        
        create_samplesheet "${batch_info}"
    fi
    
    # Start nf-core/rnaseq in the background
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting nf-core/rnaseq for ${batch_id} (in background)"

    mkdir -p "results/${batch_id}"
    touch "${batch_running_marker}"
    
    if nfcore_job_id=$(sbatch --parsable \
        --export=ALL,BATCH_ID="${batch_id}",SAMPLESHEET="${samplesheet}",BATCH_FASTQ_DIR="${batch_fastq_dir}",BATCH_DONE_MARKER="${batch_done_marker}",BATCH_RUNNING_MARKER="${batch_running_marker}",BATCH_FAILED_MARKER="${batch_failed_marker}" \
        run_nfcore_batch.sbatch); then
        nfcore_job_ids+=("${nfcore_job_id}")
        nfcore_batches+=("${batch_id}")
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitted nf-core Slurm job ${nfcore_job_id} for ${batch_id}"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Failed to submit nf-core Slurm job for ${batch_id}"
        rm -f "${batch_running_marker}"
        touch "${batch_failed_marker}"
    fi

done

# Wait for all submitted nf-core Slurm jobs to complete
if [ "${#nfcore_job_ids[@]}" -gt 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] All conversions submitted. Waiting for ${#nfcore_job_ids[@]} nf-core Slurm jobs to finish..."

    job_list=$(IFS=,; echo "${nfcore_job_ids[*]}")
    while squeue -j "${job_list}" -h -o "%i" 2>/dev/null | grep -q .; do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core Slurm jobs still running..."
        sleep 20
    done

    failed=0
    for idx in "${!nfcore_job_ids[@]}"; do
        job_id="${nfcore_job_ids[$idx]}"
        batch_id="${nfcore_batches[$idx]}"
        batch_done_marker="results/${batch_id}/.done"
        batch_running_marker="results/${batch_id}/.running"
        batch_failed_marker="results/${batch_id}/.failed"

        if [ -f "${batch_done_marker}" ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core Slurm job ${job_id} for ${batch_id} completed successfully"
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: nf-core Slurm job ${job_id} for ${batch_id} did not produce .done marker"
            [ -f "${batch_running_marker}" ] && rm -f "${batch_running_marker}"
            touch "${batch_failed_marker}"
            failed=1
        fi
    done

    if [ "${failed}" -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Some nf-core Slurm jobs failed"
        exit 1
    fi
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All batches completed successfully."
