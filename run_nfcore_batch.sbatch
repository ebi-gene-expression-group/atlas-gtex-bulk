#!/usr/bin/env bash
#SBATCH --job-name=gtex_nfcore_batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=logs/gtex_nfcore_batch_%x_%j.out
#SBATCH --error=logs/gtex_nfcore_batch_%x_%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

BATCH_ID="${BATCH_ID:?BATCH_ID is required}"
SAMPLESHEET="${SAMPLESHEET:?SAMPLESHEET is required}"
BATCH_FASTQ_DIR="${BATCH_FASTQ_DIR:?BATCH_FASTQ_DIR is required}"
BATCH_DONE_MARKER="${BATCH_DONE_MARKER:?BATCH_DONE_MARKER is required}"
BATCH_RUNNING_MARKER="${BATCH_RUNNING_MARKER:?BATCH_RUNNING_MARKER is required}"
BATCH_FAILED_MARKER="${BATCH_FAILED_MARKER:?BATCH_FAILED_MARKER is required}"

mkdir -p logs "results/${BATCH_ID}" "work_nfcore/${BATCH_ID}"

mark_failed() {
	rm -f "${BATCH_RUNNING_MARKER}"
	touch "${BATCH_FAILED_MARKER}"
}

completed=0
trap 'status=$?; if [ "${status}" -ne 0 ] || [ "${completed}" -ne 1 ]; then echo "[$(date "+%Y-%m-%d %H:%M:%S")] ERROR: nf-core/rnaseq failed for ${BATCH_ID}; preserving FASTQ and work files for debugging."; mark_failed; fi' EXIT

RIBO_INDEX=$(grep -oP '"ribo_database_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)
RIBO_MANIFEST=$(grep -oP '"ribo_database_manifest"\s*:\s*"\K[^"]+' "conf/params.json" || true)
CONTAM_INDEX=$(grep -oP '"contamination_index"\s*:\s*"\K[^"]+' "conf/params.json" || true)

NEXTFLOW_CONFIG_ARGS=()
[ -f "conf/star.config" ] && NEXTFLOW_CONFIG_ARGS+=( -c "conf/star.config" )

module load nextflow/25.04.6

PARENT_JOB_ID="${SLURM_JOB_ID:-manual}"

nextflow run rnaseq/main.nf \
	-params-file "conf/params.json" \
	"${NEXTFLOW_CONFIG_ARGS[@]}" \
	-profile singularity \
	--input "${SAMPLESHEET}" \
	--outdir "results/${BATCH_ID}" \
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
	-with-trace "${BATCH_ID}_${PARENT_JOB_ID}_trace.tsv" \
	-work-dir "work_nfcore/${BATCH_ID}" \
	-with-tower \
	-name "nf_core_gtex_rnaseq_${BATCH_ID}_${PARENT_JOB_ID}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] nf-core/rnaseq completed successfully for ${BATCH_ID}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleaning up FASTQ files and work directory..."
rm -rf "${BATCH_FASTQ_DIR}" "work_nfcore/${BATCH_ID}" ".nextflow/history" || true
rm -f "${BATCH_RUNNING_MARKER}" "${BATCH_FAILED_MARKER}"
touch "${BATCH_DONE_MARKER}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup complete for ${BATCH_ID}; marked batch as complete."

completed=1
