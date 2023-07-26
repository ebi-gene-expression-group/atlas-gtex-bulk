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
The goal of this repo in to analyse bulk GTEx V8 data (study id: `E-GTEX-8`) with a Snakemake workflow to uncompress the bams and analyse them on the fly. This can be done only by authorised users. 

1. We want to keep same tools as in the standard Atlas RNA-seq pipeline with ISL/IRAP including QC steps to flag problematic samples, with special atention to delete Fastqs and intermediate files after successful procesing.

2. Because we need to process 17350 libraries, the workflow should have constraint to enable batch processing of few `n` libraries in parallel. 

3. Input data is in BAM format

4. Output for each library should be similar to `$IRAP_SINGLE_LIB/out`.

5. Once all libraries have been processed successfully, a final aggregation rule should write final results for E-GTEX-8 in a format similar to studies here `$IRAP_SINGLE_LIB/studies`.

### Example
`snakemake -p --use-conda --conda-frontend conda --profile lsf --prioritize prepare_aggregation --keep-going --cores 4 --restart-times 5 --latency-wait 150 --config input_path=test-data atlas_gtex_root=/repo_directory_path/ private_script=gitlab_dir irap_config=homo_sapiens.conf`

For batching, we can utilise the following batch command to run few samples at a time
`snakemake -s Snakefile --cores 2 --batch final_workflow_check=n/N` where `N` is the total number of chunks, and `n=1,2, ..N`. 

For instace, if we run the workflow in `N=347` batches, 50 libraries will be processed in each batch.


### Test data
At the moment some publicly available alignment (BAM) files are available in`test-data` directory. For further analysis of iRAP/ISL pipeline more data can be downloaded following [iRAP setup data](https://github.com/nunofonseca/irap/wiki/7-Quick-Example#setup-the-data) wiki.

### Requirements
- Snakemake 7.25.3 or higher
- [LSF profile](https://github.com/Snakemake-Profiles/lsf) configuration
- Two scripts located at the config `private_script`:
  - gtex_bulk_env.sh
  - gtex_bulk_init.sh