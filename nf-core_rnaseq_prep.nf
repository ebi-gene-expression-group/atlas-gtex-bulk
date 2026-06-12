nextflow.enable.dsl=2

params.manifest   = params.manifest   ?: "bams_manifest.csv"
params.batch_size = params.batch_size ?: 100
params.batchinfo  = params.batchinfo  ?: "batch_info"

workflow {

    Channel
        .fromPath(params.manifest)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample as String
            def bam = file(row.bam as String)
            def strandedness = row.strandedness ?: "auto"

            tuple(sample, bam, strandedness)
        }
        .ifEmpty {
            error "No rows found in manifest: ${params.manifest}"
        }
        .set { bam_ch }

    /*
     Detect PE/SE per BAM.
    */
    DETECT_READ_TYPE(bam_ch)

    DETECT_READ_TYPE.out
        .map { sample, bam, strandedness, read_type_file ->
            def read_type = read_type_file.text.trim()
            tuple(sample, bam, strandedness, read_type)
        }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .flatMap { rows ->
            rows.withIndex().collect { item, idx ->
                def batch_num = Math.floor(idx / params.batch_size) + 1
                def batch_id = String.format("batch_%06d", batch_num as int)
    
                tuple(
                    batch_id,
                    item[0],             // sample
                    item[1],             // bam
                    item[2],             // strandedness
                    item[3]              // read_type
                )
            }
        }
        .groupTuple(by: 0)
        .set { batched_bams_ch }
    
    CREATE_BATCH_INFO(batched_bams_ch)
}

process DETECT_READ_TYPE {

    tag "${sample}"
    conda "envs/samtools.yml"

    input:
    tuple val(sample), path(bam), val(strandedness)

    output:
    tuple val(sample), path(bam), val(strandedness), path("read_type.txt")

    script:
    """
    set -euo pipefail

    echo "Checking BAM integrity: ${bam}"
    samtools quickcheck -v ${bam}

    echo "Detecting read type for ${sample}"

    # Count primary reads carrying the paired-end flag without using a pipe,
    # which can misclassify PE BAMs under `set -o pipefail`.
    if [ "\$(samtools view -c -f 1 -F 0x900 ${bam})" -gt 0 ]; then
        echo "pe" > read_type.txt
    else
        echo "se" > read_type.txt
    fi

    echo "Detected read type for ${sample}: \$(cat read_type.txt)"
    """
}

process CREATE_BATCH_INFO {

    tag "${batch_id}"

    publishDir "${params.batchinfo}", mode: 'copy'

    input:
    tuple val(batch_id), val(samples), val(bams), val(strandednesses), val(read_types)

    output:
    path("${batch_id}_info.csv")

    script:
    def rows = []

    for (int i = 0; i < samples.size(); i++) {
        rows << [
            samples[i],
            bams[i],
            strandednesses[i],
            read_types[i]
        ]
    }

    def csv_lines = rows.collect { row ->
        "${row[0]},${row[1]},${row[2]},${row[3]}"
    }.join("\n")

    """
    set -euo pipefail

    cat > ${batch_id}_info.csv <<'EOF'
sample,bam,strandedness,read_type
${csv_lines}
EOF

    echo "Created ${batch_id}_info.csv with \${samples.size()} samples"
    cat ${batch_id}_info.csv
    """
}
