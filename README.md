# Atlas analysis for controlled-access GTEx bulk dataset

## Justification
Controlled access (CA) data in Atlas can be processed depending on the datasets being:
- Bulk (gxa)
  - EGA/ENA - the code for EGA analysis is currently part of ISL, while the download parts are in their own repo: https://github.com/ebi-gene-expression-group/ega_downloader
  - Non-EGA, such as [GTEx](https://gtexportal.org/home/datasets)
- Single cell (scxa)
  - There’s also provision in the single-cell pipelines (currently unused) for single-cell CA data. Here’s some background on what was set up for single-cell: https://github.com/ebi-gene-expression-group/scxa-control-workflow/pull/16
  - An alternative path for ingesting data to SCXA would be via annData with our [atlas-anndata](https://github.com/ebi-gene-expression-group/atlas-anndata) tool. For instance, metadata has been extracted from annData in the GTEx portal under accession `E-ANND-2`.

## GTEx analysis
This repository currently uses a two-stage Nextflow + SLURM workflow for controlled-access GTEx bulk BAM data.

### Current workflow
1. Prepare batch metadata from BAM manifest (`nf-core_rnaseq_prep.nf`):
  - Reads `bams_manifest.csv`
  - Detects read type per BAM
  - Writes batch files in `batch_info/` as `batch_XXXXXX_info.csv`
2. Run per-batch conversion + nf-core/rnaseq (`run_nf-core_rnaseq.sh` submitted via `sbatch`):
  - Converts BAM to FASTQ per batch (parallelized within the batch)
  - Builds batch samplesheet
  - Runs `rnaseq/main.nf`
  - Cleans FASTQ + batch work directory on success
  - Writes completion marker `results/<batch_id>/.done`

### Batch prep step
Default batch size is 100 BAMs (configured in `nextflow.config`).

Run batch preparation:

```bash
nextflow run nf-core_rnaseq_prep.nf -profile singularity -resume
```

Override batch size when needed:

```bash
nextflow run nf-core_rnaseq_prep.nf -profile singularity --batch_size 50 -resume
```

### Per-batch processing step
Submit the SLURM wrapper job:

Environment setup example (HPC):

  module load nextflow
  export SAMTOOLS_IMAGE=/path/to/samtools.sif
  export TOWER_ACCESS_TOKEN=your_tower_token
  export TOWER_WORKSPACE_ID=your_tower_workspace_id

```bash
SAMTOOLS_IMAGE=/path/to/samtools.sif sbatch run_nf-core_rnaseq.sh
```

Notes:
- SAMTOOLS_IMAGE is required by run_nf-core_rnaseq.sh for BAM to FASTQ conversion.
- TOWER_ACCESS_TOKEN and TOWER_WORKSPACE_ID are required when using -with-tower.
- If your cluster uses environment modules, load Nextflow before submitting the job.

The script will:
- Process each `batch_info/batch_*.csv`
- Skip batches that already have `results/<batch_id>/.done`
- Skip FASTQ conversion when reusable inputs already exist
- Run nf-core/rnaseq for the batch
- On successful rnaseq completion:
  - remove `fastq_batches/<batch_id>`
  - remove `work_nfcore/<batch_id>`
  - create `results/<batch_id>/.done`

### Rerun behavior
- Safe reruns are supported.
- Completed batches are skipped via `.done` marker.
- Failed batches are retried on subsequent runs.

### Logging and trace
- Batch logs include timestamps.
- Nextflow trace is written per batch as `<batch_id>_trace.tsv`.
- Main SLURM log files are written to `logs/`.

### Post-processing
After all batches complete successfully, consider sanitizing MultiQC reports:

```bash
# Sanitize MultiQC HTML/JSON reports to remove absolute paths
# and temporary work directory references
find results -name "multiqc*.html" -o -name "multiqc*.json" | while read f; do
    sed -i 's|'"$(pwd)"'|.|g' "$f"
done
```

This step:
- Replaces absolute workspace paths with relative paths (`.`)
- Removes potentially sensitive file paths from reports
- Makes reports portable across systems

### Requirements
- Nextflow
- SLURM
- Singularity/Apptainer
- Samtools Singularity image (path passed via `SAMTOOLS_IMAGE`)

### Test data
Public test BAM files are available under `test-data/` for functional checks.
