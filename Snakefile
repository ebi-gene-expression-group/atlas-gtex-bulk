import os
import glob

SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.bam")

# probably need to specify here (or in config) the species or the location of the genome/annotations

rule all:
    input: expand(["out/{sample}/{sample}_1_fastqc.html","out/{sample}/{sample}_2_fastqc.html","out/{sample}/{sample}_1.fastq.gz", "out/{sample}/{sample}_2.fastq.gz"], sample=SAMPLES)

rule bam_to_fastq:
	input:
		config["input_path"]+"/{sample}.bam",
	output:
		fastq_1="out/{sample}/{sample}_1.fastq.gz",
                fastq_2="out/{sample}/{sample}_2.fastq.gz"
	conda: 
		"envs/samtools.yml"
	shell:
		"""
			samtools fastq -c 6 -1 {output.fastq_1} -2 {output.fastq_2} {input}
		"""
 
rule qc:
	input:
		fq1 = rules.bam_to_fastq.output.fastq_1,
		fq2 = rules.bam_to_fastq.output.fastq_2,
	output:
		fqc1 = "out/{sample}/{sample}_1_fastqc.html",
		fqc2 = "out/{sample}/{sample}_2_fastqc.html",
	params:
		fqc_dir = "out/{sample}/",
	conda:
		"envs/fastqc.yml"
	shell:
		"""
			fastqc {input.fq1} --outdir={params.fqc_dir}
			fastqc {input.fq2} --outdir={params.fqc_dir}
		"""


rule run_irap:
	

rule aggreagate_libraries:
    """
    Final rule to agregate all library outputs
    """


