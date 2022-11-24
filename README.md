# Atlas analysis for controlled-access GTEx bulk dataset

Controlled access (CA) data in Atlas can be processed depending on the datasets being:
- Bulk
  - EGA/ENA - the code for EGA analysis is currently part of ISL, while the download parts are in their own repo: https://github.com/ebi-gene-expression-group/ega_downloader
  - Non-EGA, such as [GTEx](https://gtexportal.org/home/datasets).
- Single cell
  - There’s also provision in the single-cell pipelines (currently unused) for single-cell CA data. Here’s some background on what was set up for single-cell: https://github.com/ebi-gene-expression-group/scxa-control-workflow/pull/16
  - An alternative path for ingesting data to SCXA would be via annData with our [atlas-anndata](https://github.com/ebi-gene-expression-group/atlas-anndata) tool. For instance, metadata has been extracted from annData in the GTEx portal under accession `E-ANND-2`.

The goal of this repo in to analyse bulk GTEx V8 bulk data (study id: `E-GTEX-8`) with a Snakemake workflow to uncompress the bams and analyse on them fly. This can be done only by authorised users. 

1. We want to keep same tools as in the standard Atlas RNA-seq pipeline with ISL/IRAP including QC steps to flag problematic samples, with special atention to delete Fastqs and intermediate files after successful procesing.

2. Because we need to process ~17300 libraries, the workflow should have constraint to enable batch processing of few libraries in parallel. 

3. Input data is in BAM format

4. Output for each library should be similar to `$IRAP_SINGLE_LIB/out`.

5. Once all libraries have been processed successfully, a final aggregation rule should write final results for E-GTEX-8 in a format similar to studies here `$IRAP_SINGLE_LIB/`.

## Example
`snakemake -p --use-conda --cores 2 --config input_path=test-data`
