# Atlas analysis for controlled-access GTEx bulk datasets

Controlled access (CA) data in Atlas can be processed depending on the datasets being:
- EGA/ENA - the code for EGA analysis is currently part of ISL, while the download parts are in their own repo: https://github.com/ebi-gene-expression-group/ega_downloader
- Non-EGA, such as GTEx.

The goal of this repo in to analyse bulk GTEx data with a Snakemake workflow to uncompress the bams and analyse on the fly. This can be done only by authorised users.

There’s also provision in the single-cell pipelines (currently unused) for single-cell CA data. Here’s some background on what was set up for single-cell: https://github.com/ebi-gene-expression-group/scxa-control-workflow/pull/16
