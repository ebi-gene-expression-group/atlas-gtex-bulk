#!/usr/bin/env bash
#SBATCH --job-name=gtex_sample_convert
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gtex_sample_convert_%j.out
#SBATCH --error=logs/gtex_sample_convert_%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

SAMPLE="${SAMPLE:?SAMPLE is required}"
BAM="${BAM:?BAM is required}"
BATCH_FASTQ_DIR="${BATCH_FASTQ_DIR:?BATCH_FASTQ_DIR is required}"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:?SAMTOOLS_IMAGE is required}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p logs

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converting ${SAMPLE} from ${BAM}"

singularity exec "${SAMTOOLS_IMAGE}" samtools collate \
    -@ "${THREADS}" \
    -u \
    -O "${BAM}" \
| singularity exec "${SAMTOOLS_IMAGE}" samtools fastq \
    -@ "${THREADS}" \
    -F 0x900 \
    -1 "${BATCH_FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
    -2 "${BATCH_FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n -

if [ $? -eq 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully converted ${SAMPLE}"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Conversion failed for ${SAMPLE}"
    exit 1
fi
