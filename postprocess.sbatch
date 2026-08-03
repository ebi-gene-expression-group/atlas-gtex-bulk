#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

set -euo pipefail
set -x 


bam_list=$1
g=$2
chrom_sizes=$3

bws=()

: "${BED_CONDA_ENV:?BED_CONDA_ENV is not set}"

env_dir="$BED_CONDA_ENV"

export PATH="${env_dir}/bin:${PATH}"

while IFS= read -r bam; do
    [[ -z "$bam" ]] && continue
    bws+=("bw/$(basename "$bam" .bam).CPM.bw")
done < "$bam_list"

N=${#bws[@]}

echo "Number of bigWigs: $N"

echo "Merging bigWigs"

bigWigMerge \
    "${bws[@]}" \
    "bw/${g}.bedGraph"

echo "Calculating mean bedGraph"

awk -v n="$N" '
    BEGIN {
        OFS = "\t"
    }
    {
        $4 = $4 / n
        print
    }
' "bw/${g}.bedGraph" > "bw/${g}.mean.bedGraph"

rm bw/${g}.bedGraph

echo "Sorting mean bedGraph"

LC_ALL=C sort \
    --parallel="$SLURM_CPUS_PER_TASK" \
    -k1,1 \
    -k2,2n \
    "bw/${g}.mean.bedGraph" > "bw/${g}.mean.expressions.bedGraph"

rm bw/${g}.mean.bedGraph

echo "Creating mean bigWig"

bedGraphToBigWig \
    "bw/${g}.mean.expressions.bedGraph" \
    "$chrom_sizes" \
    "bw/${g}.mean.CPM.bw"


echo "Creating D4"

d4tools create \
    "bw/${g}.mean.CPM.bw" \
    "bw/${g}.mean.CPM.d4"

touch "bw/${g}.done"

rm -rf ${bws[@]}

echo "Completed group: $g"
