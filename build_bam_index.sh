#!/usr/bin/env bash
set -euo pipefail

list_dir="assays"
search_root="results"
index_file="bam_index.tsv"

: > "$index_file"

for bam in "$search_root"/batch_*/star_salmon/*.bam; do
    [[ -f "$bam" ]] || continue

    filename=${bam##*/}
    printf '%s\t%s\n' "$filename" "$bam"
done > "$index_file"

printf 'Indexed %s BAM files\n' "$(wc -l < "$index_file")" >&2
