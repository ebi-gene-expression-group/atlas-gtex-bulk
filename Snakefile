import os
import glob

# here we need wildcard contraint based on input file name pattern or it picks up bam file from sub-dirs
SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.bam")

# probably need to specify here (or in config) the species or the location of the genome/annotations

rule all:
    input: expand(["out/{sample}/{sample}.fastq.val", "out/{sample}.txt"], sample=SAMPLES)


rule check_bam:
    input:
        config["input_path"]+ "/{sample}.bam"
    output:
        bam_check = "out/{sample}/{sample}_1.bam_checked",
    conda:
        "envs/samtools.yml"
    shell:
        """
        samtools quickcheck {input}
	if [ $? -ne 0 ]; then
		echo "{input} is not valid BAM"
	else
		touch {output.bam_check}
        fi
        """

rule bam_to_fastq:
    input:
        bam = config["input_path"]+"/{sample}.bam",
        check_bam = rules.check_bam.output.bam_check
    output:
        fastq = "out/{sample}/{sample}.fq.gz",
    conda: 
        "envs/samtools.yml"
    shell:
        """
        samtools fastq -0 /dev/null {input.bam} > {output.fastq}
        """

# here if else statement could be modified to run iRAP/ISL for PE and SE data
checkpoint validating_fastq:
    input:
        fastq = rules.bam_to_fastq.output.fastq,
    output:
        val_fastq = "out/{sample}/{sample}.fastq.val",
    params:
        fastq = "out/{sample}/{sample}.fastq.gz",
    conda:
        "envs/fastq_utils.yml"
    shell:
        """
        fastq_info {input.fastq}  
        if [ $? -ne 0 ]; then
            #rm -rf {input.fastq}
            echo "ERROR: Failed fastq validation {input.fastq}"
            mv {input.fastq} {params.fastq}
        else
            echo "validation successful"
        fi

        touch {output.val_fastq}
        """

def aggregate_fastq(wildcards):
    checkpoint_output = checkpoints.validating_fastq.get(**wildcards).output[0]
    print(checkoint_output)
    return expand("out/{sample}/{fq}_fastqc.html",
        sample=wildcards.sample,
        fq=glob_wildcards(os.path.join(checkpoint_output, "{fq}.fastq.gz")).fq)

rule qc:
    input:
        check = rules.validating_fastq.output.val_fastq
    output:
        fqc = "out/{sample}/{fq}_1_fastqc.html",
    params:
        fq = "out/{sample}/{fq}.fastq.gz",
        fqc_dir = "out/{sample}"
    conda:
        "envs/fastqc.yml"
    shell:
        """
        fastqc {params.fq} --outdir={params.fqc_dir}
        """

rule run_irap:
        input:
                fastq=rules.bam_to_fastq.output.fastq,
        output: "out/{sample}.txt"
	conda: "envs/isl.yaml"
        params:
                private_script=config["private_script"],
                conf=config["irap_config"],
                strand="both",
                irapMem=4096000000,
                irapDataOption="",
		filename="{sample}"
	resources: mem_mb=10000	
	shell:
                """
                source {params.private_script}/gtex_bulk_env.sh
                source {params.private_script}/gtex_bulk_init.sh 
                source {params.private_script}/irap.sh
		cp {params.private_script}/gtex_bulk_env.sh $IRAP_SINGLE_LIB
		cat {input.fastq} | grep '^@.*/1$' -A 3 --no-group-separator > {params.filename}_1.fastq
		cat {input.fastq} | grep '^@.*/2$' -A 3 --no-group-separator > {params.filename}_2.fastq
                fastq_info {params.filename}_1.fastq {params.filename}_2.fastq
	        if [ $? -ne 0 ]; then
	        #fastq is SE
	        # iRAP SE command here
	        echo "SE "
	        else
	        # fastq is PE
		# iRAP PE commands here
		cmd="irap_single_lib -A -f -o irap_single_lib/{params.filename} -1 {params.filename}_1.fastq -2 {params.filename}_2.fastq -c {params.conf} -s {params.strand} -m {params.irapMem} -t 5 -C {params.irapDataOption}"
		echo "PE"
                eval $cmd
                fi
		touch {output}
"""


rule isl_db_update:
    """
    In manual processing we need to update the ISL LIBRARIES table
    DB cannot be accessed by user other than fg_atlas, 
    so leave a csv file (or lock files) of updates somewhere that ISL can access
    and use to update the dbs during usual processing.
    This will need additional logic in repo `isl`.
    """

rule aggreagate_libraries:
    """
    Final rule to agregate all library outputs, and store them in a single folder
    $IRAP_SINGLE_LIB/studies/E-GTEX-8/homo_sapiens/
    """

rule final_check:
    """
    should look for .complete file
    """
