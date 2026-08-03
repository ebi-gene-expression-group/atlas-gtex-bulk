# Atlas analysis for controlled-access GTEx bulk dataset

## Justification
Controlled-access (CA) data in Atlas can be processed depending on the dataset type:

- Bulk (`gxa`)
  - EGA/ENA: EGA analysis code is currently part of ISL, while the download code is in https://github.com/ebi-gene-expression-group/ega_downloader
  - Non-EGA controlled-access data, such as [GTEx](https://gtexportal.org/home/datasets)
- Single cell (`scxa`)
  - There is currently unused provision in the single-cell pipelines for CA data. Background: https://github.com/ebi-gene-expression-group/scxa-control-workflow/pull/16
  - An alternative path for SCXA ingest is via AnnData using [atlas-anndata](https://github.com/ebi-gene-expression-group/atlas-anndata). For example, GTEx portal metadata has been extracted from AnnData under accession `E-ANND-2`.

## Repository Overview
This repository runs controlled-access GTEx bulk BAMs through an nf-core/rnaseq-style workflow on SLURM.

Main files:

- `bams_manifest.csv`: input manifest with `sample,bam,strandedness` columns.
- `nf-core_rnaseq_prep.nf`: validates BAMs, detects paired-end/single-end read type, and writes per-batch CSVs.
- `submit_nfcore_rnaseq_prep.sbatch`: SLURM wrapper for batch metadata preparation.
- `nextflow.config`: SLURM/Singularity config for the prep workflow and selected nf-core process resources.
- `run_nf-core_rnaseq.sh`: top-level SLURM wrapper for all prepared batches.
- `convert_sample_fastq.sbatch`: converts one BAM to paired FASTQ files.
- `convert_batch_fastq.sbatch`: optional standalone per-batch FASTQ conversion helper.
- `run_nfcore_batch.sbatch`: runs `rnaseq/main.nf` for one prepared batch.
- `conf/params.json`: nf-core/rnaseq reference and contamination/ribo database parameters.
- `conf/star.config`: STAR argument override used automatically if present.
- `rnaseq/`: nf-core/rnaseq pipeline submodule.
- `merge_existing_star_salmon_matrices.py`: merges selected STAR-Salmon matrices from completed batch result directories.
- `merge_nfcore_rnaseq_batches.py`: creates a broader nf-core-like combined output from multiple batch result directories.
- `create_star_salmon_merged_metrices.sh`: SLURM wrapper with a hard-coded full-batch merge command; edit the batch list before use.
- `download/`: older GTEx download helper scripts and notes.
- `test-data/`: public mini BAMs for functional checks.

If cloning fresh, initialise the pipeline submodule:

```bash
git submodule update --init --recursive
```

## Requirements
- Nextflow
- SLURM with `sbatch` and `squeue`
- Singularity/Apptainer
- A Samtools Singularity image, passed as `SAMTOOLS_IMAGE`
- Tower/Seqera credentials if using the default `run_nfcore_batch.sbatch`, which runs nf-core/rnaseq with `-with-tower`
- Python with pandas for merge scripts; `envs/nfcore-merge.yaml` provides a conda environment

The batch nf-core wrapper currently loads `nextflow/25.04.6` inside `run_nfcore_batch.sbatch`. Adjust that module line if your cluster uses a different module name.

## Configuration
### BAM Manifest
`bams_manifest.csv` must contain:

```csv
sample,bam,strandedness
sample_1,/path/to/sample_1.bam,auto
```

`strandedness` is passed through to the nf-core samplesheet. Use values accepted by nf-core/rnaseq, for example `auto`, `forward`, `reverse`, or `unstranded`.

### nf-core Parameters
Update `conf/params.json` before running on real data. It contains paths for:

- Genome FASTA/GTF
- STAR and Salmon indexes
- Contamination Kraken2/Bracken database
- SortMeRNA ribosomal RNA manifest and index

Use real absolute paths in `conf/params.json`. The shell wrapper extracts some values from this JSON and passes them as command-line arguments, so placeholder strings such as `${BULK_REFERENCES_DIR}` will not be expanded by the shell at that point.

`conf/star.config` is optional, but is included automatically by `run_nfcore_batch.sbatch` when present.

## Workflow
### 1. Prepare Batch Metadata
Default batch size is 100 BAMs, configured in `nextflow.config`.

Run directly:

```bash
nextflow run nf-core_rnaseq_prep.nf -c nextflow.config -profile singularity -work-dir work_nfprep -resume
```

Or submit the SLURM wrapper:

```bash
sbatch submit_nfcore_rnaseq_prep.sbatch
```

Useful overrides:

```bash
BATCH_SIZE=50 MANIFEST=bams_manifest.csv sbatch submit_nfcore_rnaseq_prep.sbatch
```

This step writes:

```text
batch_info/
  batch_000001_info.csv
  batch_000002_info.csv
```

Each batch info file contains:

```csv
sample,bam,strandedness,read_type
```

### 2. Run All Prepared Batches
Submit the top-level wrapper:

```bash
export SAMTOOLS_IMAGE=/path/to/samtools.sif
export TOWER_ACCESS_TOKEN=your_tower_token
export TOWER_WORKSPACE_ID=your_tower_workspace_id

sbatch run_nf-core_rnaseq.sh
```

For each `batch_info/batch_*.csv`, the wrapper:

- Skips batches with `results/<batch_id>/.done`, `.running`, or `.failed`.
- Converts each BAM with `convert_sample_fastq.sbatch`, unless the expected FASTQs and samplesheet already exist.
- Writes `fastq_batches/<batch_id>/samplesheet.csv`.
- Submits `run_nfcore_batch.sbatch` for the batch.
- Waits for submitted nf-core jobs to finish.
- Requires `results/<batch_id>/.done` to consider a batch successful.

`run_nfcore_batch.sbatch` then:

- Runs `rnaseq/main.nf` with `conf/params.json`.
- Adds `conf/star.config` via `-c` if present.
- Writes results to `results/<batch_id>`.
- Writes work files to `work_nfcore/<batch_id>`.
- Writes a trace file named `<batch_id>_<parent_slurm_job_id>_trace.tsv`.
- Removes the batch FASTQs and `work_nfcore/<batch_id>` after successful completion.
- Creates `results/<batch_id>/.done` on success.
- Creates `results/<batch_id>/.failed` on failure and preserves FASTQs/work files for debugging.

### Optional Manual Per-Batch Run
`convert_batch_fastq.sbatch` can convert an entire prepared batch in one SLURM job:

```bash
BATCH_INFO=batch_info/batch_000001_info.csv \
SAMTOOLS_IMAGE=/path/to/samtools.sif \
sbatch convert_batch_fastq.sbatch
```

Then submit the nf-core batch job manually:

```bash
BATCH_ID=batch_000001 \
SAMPLESHEET=fastq_batches/batch_000001/samplesheet.csv \
BATCH_FASTQ_DIR=fastq_batches/batch_000001 \
BATCH_DONE_MARKER=results/batch_000001/.done \
BATCH_RUNNING_MARKER=results/batch_000001/.running \
BATCH_FAILED_MARKER=results/batch_000001/.failed \
sbatch run_nfcore_batch.sbatch
```

## Rerun Behavior
- Completed batches are skipped with `.done`.
- Running batches are skipped with `.running`.
- Failed batches are skipped with `.failed`.
- To retry a failed batch, inspect the preserved files, then remove `results/<batch_id>/.failed` and rerun `run_nf-core_rnaseq.sh`.
- Existing FASTQs are reused when both the samplesheet and expected FASTQ files are present.
- Successful batches clean up `fastq_batches/<batch_id>` and `work_nfcore/<batch_id>`.

## Logs and Outputs
- SLURM logs are written to `logs/`.
- Batch FASTQs and generated samplesheets are written under `fastq_batches/` until successful cleanup.
- nf-core/rnaseq batch outputs are written to `results/<batch_id>/`.
- nf-core work directories are written under `work_nfcore/<batch_id>/`.
- Prep workflow work directories are written under `work_nfprep/`.
- Nextflow trace files are written as `<batch_id>_<slurm_job_id>_trace.tsv`.

## Post-Processing
After all batches complete successfully, consider sanitizing MultiQC reports to remove absolute workspace paths:

```bash
find results \( -name "multiqc*.html" -o -name "multiqc*.json" \) -type f | while read -r f; do
    sed -i 's|'"$(pwd)"'|.|g' "$f"
done
```

This replaces absolute workspace paths with `.` and makes reports more portable.

## Aggregation
Create the merge environment if needed:

```bash
conda env create -f envs/nfcore-merge.yaml
conda activate nfcore-merge
```

To merge the four existing STAR-Salmon matrix files from completed batch directories:

```bash
python merge_existing_star_salmon_matrices.py \
  -o batch_1-2 \
  results/batch_000001 \
  results/batch_000002
```

This writes:

```text
batch_1-2/
└── star_salmon
    ├── salmon.merged.gene_counts.tsv
    ├── salmon.merged.gene_tpm.tsv
    ├── salmon.merged.transcript_counts.tsv
    └── salmon.merged.transcript_tpm.tsv
```

Use `--relaxed` only when batches do not have identical feature rows/order and an outer join is acceptable.

For a broader nf-core-like combined directory with merged known matrices, linked/copied sample directories, metadata, and a manifest:

```bash
python merge_nfcore_rnaseq_batches.py \
  -o batch_all \
  results/batch_000001 \
  results/batch_000002
```

Add `--copy` to copy sample directories instead of symlinking them. Add `--relaxed-features` only when feature IDs/order differ and an outer join is acceptable.

`create_star_salmon_merged_metrices.sh` is a convenience SLURM script for merging many batch directories into `batch_all`; review and edit its hard-coded batch list before submitting it.

## Expression coverage generation
Create the bed and bw environments if needed:

```bash
conda env create -f envs/bed.yaml
conda env create -f envs/bw.yaml
```

To create expression coverage 
1. First we create index for alignment files (.bam) using build_bam_index.sh which will create bam_index.tsv containing alignment file and path.
2. Create file per assay containing samples.
    ```
      xml_file=E-GTEX-8-configuration.xml
      output_dir=assays
      awk -v out="$output_dir" '
      /<assay_group[[:space:]][^>]*id="/ {
          group_id = $0
          sub(/^.*id="/, "", group_id)
          sub(/".*$/, "", group_id)
      }
      /<assay>/ {
          assay = $0
          sub(/^.*<assay>[[:space:]]*/, "", assay)
          sub(/[[:space:]]*<\/assay>.*$/, "", assay)
          if (group_id != "" && assay != "")
              print assay >> out "/" group_id ".txt"
      }
      /<\/assay_group>/ {
          if (group_id != "")
              close(out "/" group_id ".txt")
          group_id = ""
      }
      ' "$xml_file"
    ```
3. Create symlink of fa.sizes file.  `ln -sf <path>/references/homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.toplevel.fa.sizes Homo_sapiens.GRCh38.dna.toplevel.fa.sizes`
4. Export environment variable `BED_CONDA_ENV` and `BW_CONDA_ENV` pointing to conda env path
5. Generate bw and d4 per assay by running `bas submit_groups.sh`

  ```bash
    python merge_existing_star_salmon_matrices.py \
      -o batch_1-2 \
      results/batch_000001 \
      results/batch_000002
  ```

This writes:


## Functional Check
The repository includes public mini BAMs under `test-data/`, and the checked-in `bams_manifest.csv` points at them. A small local/HPC smoke test is:

```bash
nextflow run nf-core_rnaseq_prep.nf -c nextflow.config -profile singularity --batch_size 2 -work-dir work_nfprep -resume
```

This should create small `batch_info/batch_XXXXXX_info.csv` files. Full nf-core execution still requires valid references in `conf/params.json`, a Samtools image, and the `rnaseq/` submodule.
