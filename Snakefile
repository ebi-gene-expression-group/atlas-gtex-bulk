import os
import glob

SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.bam")

# probably need to specify here (or in config) the species or the location of the genome/annotations

rule all:
    input: expand(["out/{sample}/{sample}.fastq.val"], sample=SAMPLES)


rule check_bam:
    input:
        config["input_path"]+ "/{sample}.bam"
    output:
        bam_check = temp("out/{sample}/{sample}_1.bam_checked"),
        test_out = temp("out/{sample}/{sample}_1_check.txt"),
        test_log = temp("out/{sample}/{sample}_1_check.log")
    conda:
        "envs/rsamtools.yml"
    shell:
        """
        Rscript scripts/testPairedEndBam.R {input} > {output.test_out} 2> {output.test_log}
        if cat {output.test_out} | grep -q "TRUE"; then
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
        else
            echo "validation successful"
            mv {input.fastq} {params.fastq}
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
    """
    Run IRAP analysis. Need conda isl env and singularity
    """
    #input:
    #    fastq_1 = rules.qc.output.fqc,
    #params:
        #irap_container=config['irap_container']
    #conda:
    	#"envs/isl.yml"   ### change to submodule env once ISL is linked as submodule
    #output:


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
