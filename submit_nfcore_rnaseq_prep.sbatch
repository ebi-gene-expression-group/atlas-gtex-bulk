#!/usr/bin/env bash
#SBATCH --job-name=nfcore_prep
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/nfcore_prep.%j.out
#SBATCH --error=logs/nfcore_prep.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your.email@example.com

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

mkdir -p logs .nextflow work_nfprep

PIPELINE_SCRIPT="${PIPELINE_SCRIPT:-nf-core_rnaseq_prep.nf}"
NF_CONFIG="${NF_CONFIG:-nextflow.config}"
PROFILE="${PROFILE:-singularity}"
MANIFEST="${MANIFEST:-bams_manifest.csv}"
BATCH_SIZE="${BATCH_SIZE:-100}"
FASTQ_OUTDIR="${FASTQ_OUTDIR:-fastq_batches}"
SAMPLESHEETS="${SAMPLESHEETS:-samplesheets}"
USE_RESUME="${USE_RESUME:-true}"

cmd=(
  nextflow run "${PIPELINE_SCRIPT}"
  -c "${NF_CONFIG}"
  -profile "${PROFILE}"
  --manifest "${MANIFEST}"
  --batch_size "${BATCH_SIZE}"
  --fastq_outdir "${FASTQ_OUTDIR}"
  --samplesheets "${SAMPLESHEETS}"
  -work-dir "work_nfprep"
)

if [[ "${USE_RESUME}" == "true" ]]; then
  cmd+=( -resume )
fi

"${cmd[@]}"
