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

declare -a nfcore_pids=()
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
    samplesheet="${batch_fastq_dir}/samplesheet.csv"
    
    echo "=========================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${batch_id}"
    echo "=========================================="
    
    if [ -f "${batch_done_marker}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch ${batch_id} already completed; skipping."
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
    
    RIBO_INDEX=$(grep -oP '"ribo_database_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)
    RIBO_MANIFEST=$(grep -oP '"ribo_database_manifest"\s*:\s*"\K[^"]+' "conf/params.json" || true)
    CONTAM_INDEX=$(grep -oP '"contamination_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)

    NEXTFLOW_CONFIG_ARGS=()
    [ -f "conf/rnaseq.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/rnaseq.config" )
    [ -f "conf/star.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/star.config" )

    module load nextflow/25.04.6
    
    (
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
            -with-trace "${batch_id}_${SLURM_JOB_ID}_trace.tsv" \
            -work-dir "work_nfcore/${batch_id}" \
            -with-tower \
            -name "nf_core_gtex_rnaseq_${batch_id}_${SLURM_JOB_ID}"
        
        if [ $? -eq 0 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core/rnaseq completed successfully for ${batch_id}"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleaning up FASTQ files and work directory..."
            rm -rf "${batch_fastq_dir}" "work_nfcore/${batch_id}" ".nextflow/history" || true
            mkdir -p "results/${batch_id}"
            touch "${batch_done_marker}"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup complete for ${batch_id}; marked batch as complete."
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: nf-core/rnaseq failed for ${batch_id}; preserving FASTQ and work files for debugging."
            exit 1
        fi
    ) &
    
    nfcore_pids+=("$!")
    nfcore_batches+=("${batch_id}")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core job for ${batch_id} running in background (PID: $!)"

done

# Wait for all background nf-core jobs to complete
if [ "${#nfcore_pids[@]}" -gt 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] All conversions submitted. Waiting for ${#nfcore_pids[@]} nf-core jobs to finish..."
    
    failed=0
    for idx in "${!nfcore_pids[@]}"; do
        pid="${nfcore_pids[$idx]}"
        batch_id="${nfcore_batches[$idx]}"
        
        if ! wait "${pid}"; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: nf-core job for ${batch_id} (PID: ${pid}) failed"
            failed=1
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core job for ${batch_id} completed successfully"
        fi
    done
    
    if [ "${failed}" -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Some nf-core jobs failed"
        exit 1
    fi
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All batches completed successfully."
