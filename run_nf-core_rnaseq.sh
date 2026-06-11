#!/usr/bin/env bash
set -euo pipefail

NFCORE_VERSION="3.18.0"
PROFILE="singularity,slurm"
GENOME="GRCh38"
ALIGNER="star_salmon"

mkdir -p logs results work_nfcore

for ss in samplesheets/samplesheet_batch_*.csv; do

    batch=$(basename "$ss" .csv)
    batch=${batch#samplesheet_}

    echo "Submitting nf-core/rnaseq for ${batch}"

    sbatch \
      --job-name="rnaseq_${batch}" \
      --cpus-per-task=1 \
      --mem=4G \
      --time=10-00:00:00 \
      --output="logs/${batch}.%j.out" \
      --error="logs/${batch}.%j.err" \
      <<EOF
#!/usr/bin/env bash
set -euo pipefail

nextflow run nf-core/rnaseq \\
  -r ${NFCORE_VERSION} \\
  -profile ${PROFILE} \\
  --input "${ss}" \\
  --outdir "results/${batch}" \\
  --genome ${GENOME} \\
  --aligner ${ALIGNER} \\
  -work-dir "work_nfcore/${batch}" \\
  -resume

EOF

done
