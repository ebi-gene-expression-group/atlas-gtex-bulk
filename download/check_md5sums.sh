topdir=$(pwd)
file_dir='data/rnaseq_aligned_reads'
manifest='manifests/rnaseq_aligned_reads.json'

# First derive what the cheksums should be

if [ ! -e target_checksums.txt ]; then
    for row in $(cat $manifest | jq -r '.[] | @base64'); do
        _jq() {
         echo ${row} | base64 --decode | jq -r ${1}
        }
        echo "$(_jq '.md5sum')  $(_jq '.file_name')"
    done > target_checksums.txt
fi

# Split the checksums file to parallelise the checks

checksum_parts_dir=$(pwd)/target_checksum_parts
mkdir -p $checksum_parts_dir
rm -f ${checksum_parts_dir}/*
rm ${topdir}/md5_fails.txt
split --lines=250 target_checksums.txt ${checksum_parts_dir}/

# Change to the directory with the files and bsub the check jobs
pushd $file_dir > /dev/null 2>&1
ls $checksum_parts_dir | while read -r l; do
    bsub "md5sum -c $checksum_parts_dir/$l 2>/dev/null | grep FAILED >> ${topdir}/md5_fails.txt"
    sleep 1
done
popd > /dev/null 2>&1

