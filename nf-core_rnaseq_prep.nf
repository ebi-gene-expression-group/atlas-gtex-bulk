nextflow.enable.dsl=2

params.manifest     = params.manifest     ?: "bams_manifest.csv"
params.batch_size   = params.batch_size   ?: 100
params.fastq_outdir = params.fastq_outdir ?: "fastq_batches"
params.samplesheets = params.samplesheets ?: "samplesheets"

workflow {

    def fastq_outdir_abs = file(params.fastq_outdir).toAbsolutePath().toString()

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
                    item[3],             // read_type
                    fastq_outdir_abs
                )
            }
        }
        .set { batched_bams_ch }
    
    batched_bams_ch
        .branch {
            pe: it[4] == "pe"
            se: it[4] == "se"
        }
        .set { readtype_ch }
    
    BAM_TO_FASTQ_PE(readtype_ch.pe)
    BAM_TO_FASTQ_SE(readtype_ch.se)
    
    BAM_TO_FASTQ_PE.out
        .map { batch_id, sample, fastq_1, fastq_2, strandedness, read_type, r1_file, r2_file ->
            tuple(batch_id, sample, fastq_1, fastq_2, strandedness, read_type)
        }
        .mix(
            BAM_TO_FASTQ_SE.out.map { batch_id, sample, fastq_1, fastq_2, strandedness, read_type, fastq_file ->
                tuple(batch_id, sample, fastq_1, fastq_2, strandedness, read_type)
            }
        )
        .groupTuple(by: 0)
        .set { grouped_fastqs_ch }
    
    MAKE_SAMPLESHEET(grouped_fastqs_ch)
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




process BAM_TO_FASTQ_PE {

    tag "${batch_id}:${sample}:PE"
    conda "envs/samtools.yml"

    publishDir { "${fastq_outdir}/${batch_id}" }, mode: 'copy'

    input:
    tuple val(batch_id), val(sample), path(bam), val(strandedness), val(read_type), val(fastq_outdir)

    output:
    tuple(
        val(batch_id),
        val(sample),
        val("${fastq_outdir}/${batch_id}/${sample}_R1.fastq.gz"),
        val("${fastq_outdir}/${batch_id}/${sample}_R2.fastq.gz"),
        val(strandedness),
        val("pe"),
        path("${sample}_R1.fastq.gz"),
        path("${sample}_R2.fastq.gz")
    )

    script:
    """
    set -euo pipefail

    echo "Converting paired-end BAM to FASTQ: ${sample}"

    samtools quickcheck -v ${bam}

    samtools collate \\
        -@ ${task.cpus} \\
        -u \\
        -O ${bam} \\
    | samtools fastq \\
        -@ ${task.cpus} \\
        -F 0x900 \\
        -1 ${sample}_R1.fastq.gz \\
        -2 ${sample}_R2.fastq.gz \\
        -0 /dev/null \\
        -s /dev/null \\
        -n -

    gzip -t ${sample}_R1.fastq.gz
    gzip -t ${sample}_R2.fastq.gz

    R1_LINES=\$(gzip -cd ${sample}_R1.fastq.gz | wc -l)
    R2_LINES=\$(gzip -cd ${sample}_R2.fastq.gz | wc -l)

    if [ \$((R1_LINES % 4)) -ne 0 ]; then
        echo "ERROR: R1 FASTQ line count is not divisible by 4"
        exit 1
    fi

    if [ \$((R2_LINES % 4)) -ne 0 ]; then
        echo "ERROR: R2 FASTQ line count is not divisible by 4"
        exit 1
    fi

    R1_READS=\$((R1_LINES / 4))
    R2_READS=\$((R2_LINES / 4))

    if [ "\$R1_READS" -ne "\$R2_READS" ]; then
        echo "ERROR: PE read counts differ for ${sample}"
        echo "R1 reads: \$R1_READS"
        echo "R2 reads: \$R2_READS"
        exit 1
    fi

    echo "PE FASTQ validation passed for ${sample}"
    """
}

process BAM_TO_FASTQ_SE {

    tag "${batch_id}:${sample}:SE"
    conda "envs/samtools.yml"

    publishDir { "${fastq_outdir}/${batch_id}" }, mode: 'copy'

    input:
    tuple val(batch_id), val(sample), path(bam), val(strandedness), val(read_type), val(fastq_outdir)

    output:
    tuple(
        val(batch_id),
        val(sample),
        val("${fastq_outdir}/${batch_id}/${sample}.fastq.gz"),
        val(""),
        val(strandedness),
        val("se"),
        path("${sample}.fastq.gz")
    )

    script:
    """
    set -euo pipefail

    echo "Converting single-end BAM to FASTQ: ${sample}"

    samtools quickcheck -v ${bam}

    samtools fastq \\
        -@ ${task.cpus} \\
        -F 0x900 \\
        -n \\
        -o ${sample}.fastq.gz \\
        ${bam}

    gzip -t ${sample}.fastq.gz

    SE_LINES=\$(gzip -cd ${sample}.fastq.gz | wc -l)

    if [ \$((SE_LINES % 4)) -ne 0 ]; then
        echo "ERROR: SE FASTQ line count is not divisible by 4"
        exit 1
    fi

    SE_READS=\$((SE_LINES / 4))

    echo "SE FASTQ validation passed for ${sample}"
    echo "Reads: \$SE_READS"
    """
}

process MAKE_SAMPLESHEET {

    tag "${batch_id}"

    publishDir "${params.samplesheets}", mode: 'copy'

    input:
    tuple val(batch_id), val(samples), val(fastq_1s), val(fastq_2s), val(strandednesses), val(read_types)

    output:
    path("samplesheet_${batch_id}.csv")

    script:
    def rows = []

    for (int i = 0; i < samples.size(); i++) {
        rows << [
            samples[i],
            fastq_1s[i],
            fastq_2s[i],
            strandednesses[i]
        ]
    }

    def csv_lines = rows.collect { row ->
        "${row[0]},${row[1]},${row[2]},${row[3]}"
    }.join("\n")

    """
    set -euo pipefail

    cat > samplesheet_${batch_id}.csv <<'EOF'
sample,fastq_1,fastq_2,strandedness
${csv_lines}
EOF

    echo "Created samplesheet_${batch_id}.csv"
    cat samplesheet_${batch_id}.csv
    """
}
