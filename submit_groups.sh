#!/usr/bin/env bash
set -euo pipefail
list_dir="assays"
index_file="bam_index.tsv"
chrom_sizes="Homo_sapiens.GRCh38.dna.toplevel.fa.sizes"
postprocess_script="postprocess.sbatch"


: "${BW_CONDA_ENV:?BW_CONDA_ENV is not set}"


BAM_ENV="${BW_CONDA_ENV}"

mkdir -p logs

for list_file in "$list_dir"/*.txt; do
    [[ -f "$list_file" ]] || continue

    group_id=$(basename "$list_file" .txt)
    g="$group_id"
    bam_list="${g}.bams.txt"

    if [[ -f "bw/${g}.done" ]]; then
        echo "Skipping $group_id: ${g}.done exists"
        continue
    fi

    echo "Preparing group: $group_id"

    : > "$bam_list"

    bam_jobs=()
    declare -A seen_bams=()

    while IFS= read -r prefix || [[ -n "$prefix" ]]; do
        prefix=${prefix%$'\r'}
        [[ -z "$prefix" ]] && continue

        mapfile -t matches < <(
            awk -F '\t' -v prefix="$prefix" '
                index($1, prefix) == 1 {
                    print $2
                }
            ' "$index_file"
        )

        if (( ${#matches[@]} == 0 )); then
            echo "BAM NOT FOUND: $prefix" >&2
            continue
        fi

        for bam in "${matches[@]}"; do
            [[ -z "$bam" ]] && continue

            # Do not submit the same BAM twice.
            [[ -n ${seen_bams["$bam"]+x} ]] && continue
            seen_bams["$bam"]=1

            printf '%s\n' "$bam" >> "$bam_list"

            bw="$(basename "$bam" .bam).CPM.bw"

            #
            # Submit one BAM → bigWig job.
            #

	    bamCoverage_bin="${BAM_ENV}/bin/bamCoverage"
            bam_job=$(
                sbatch \
                    --parsable \
                    --job-name="${g}_bam2bw" \
                    --cpus-per-task=4 \
                    --mem=8G \
                    --time=12:00:00 \
                    --output="logs/${g}_bam2bw_%j.log" \
                    --wrap="$bamCoverage_bin \
                        --bam '$bam' \
                        --normalizeUsing CPM \
                        --binSize 1 \
                        --numberOfProcessors 4 \
                        --outFileName 'bw/$bw'"
            )

            bam_job=${bam_job%%;*}
            bam_jobs+=("$bam_job")

            echo "  Submitted BAM job $bam_job: $bam"
        done

    done < "$list_file"

    if (( ${#bam_jobs[@]} == 0 )); then
        echo "No BAM files found for $group_id"
        continue
    fi

    #
    # Make dependency string:
    # job1:job2:job3
    #

    dependency=$(IFS=:; echo "${bam_jobs[*]}")

    #
    # Submit external post-processing script.
    # It waits for all BAM jobs.
    #

    post_job=$(
        sbatch \
            --parsable \
            --dependency="afterok:${dependency}" \
            --job-name="${g}_postprocess" \
            --output="logs/${g}_postprocess_%j.log" \
            "$postprocess_script" \
            "$bam_list" \
            "$g" \
            "$chrom_sizes"
    )

    post_job=${post_job%%;*}

    echo "Submitted group: $group_id"
    echo "  BAM jobs: ${bam_jobs[*]}"
    echo "  Post-processing job: $post_job"
done
