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
     Assign batch IDs: batch_000001, batch_000002, ...
    */
    bam_ch
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
                    fastq_outdir_abs     // final FASTQ output root
                )
            }
        }
        .set { batched_bams_ch }

    /*
     Convert each BAM to paired FASTQ.
    */
    BAM_TO_FASTQ(batched_bams_ch)

    /*
     Drop the physical path objects before making samplesheet.
     We keep only the final published FASTQ paths.
    */
    BAM_TO_FASTQ.out
        .map { batch_id, sample, r1_final, r2_final, strandedness, r1_file, r2_file ->
            tuple(batch_id, sample, r1_final, r2_final, strandedness)
        }
        .groupTuple(by: 0)
        .set { grouped_fastqs_ch }

    /*
     Make one nf-core/rnaseq samplesheet per batch.
    */
    MAKE_SAMPLESHEET(grouped_fastqs_ch)
}


process BAM_TO_FASTQ {

    tag "${batch_id}:${sample}"

    publishDir { "${fastq_outdir}/${batch_id}" }, mode: 'copy'

    input:
    tuple val(batch_id), val(sample), path(bam), val(strandedness), val(fastq_outdir)

    output:
    tuple(
        val(batch_id),
        val(sample),
        val("${fastq_outdir}/${batch_id}/${sample}_R1.fastq.gz"),
        val("${fastq_outdir}/${batch_id}/${sample}_R2.fastq.gz"),
        val(strandedness),
        path("${sample}_R1.fastq.gz"),
        path("${sample}_R2.fastq.gz")
    )

    script:
    """
    set -euo pipefail

    echo "Checking BAM: ${bam}"
    samtools quickcheck -v ${bam}

    echo "Converting BAM to FASTQ for sample: ${sample}"

    # Exclude secondary and supplementary alignments:
    # 0x100 = secondary
    # 0x800 = supplementary
    # 0x900 = both
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

    echo "Validating FASTQ gzip integrity"
    gzip -t ${sample}_R1.fastq.gz
    gzip -t ${sample}_R2.fastq.gz

    echo "Counting FASTQ records"

    R1_LINES=\$(zcat ${sample}_R1.fastq.gz | wc -l)
    R2_LINES=\$(zcat ${sample}_R2.fastq.gz | wc -l)

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
        echo "ERROR: R1 and R2 read counts differ"
        echo "R1 reads: \$R1_READS"
        echo "R2 reads: \$R2_READS"
        exit 1
    fi

    echo "FASTQ validation passed for ${sample}"
    echo "Reads: \$R1_READS"
    """
}


process MAKE_SAMPLESHEET {

    tag "${batch_id}"

    publishDir "${params.samplesheets}", mode: 'copy'

    input:
    tuple val(batch_id), val(samples), val(r1_finals), val(r2_finals), val(strandednesses)

    output:
    path("samplesheet_${batch_id}.csv")

    script:
    def rows = []

    for (int i = 0; i < samples.size(); i++) {
        rows << [
            samples[i],
            r1_finals[i],
            r2_finals[i],
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
