import os
import glob

SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.bam")


print(SAMPLES)

rule all:
    input: expand(["out/{sample}/{sample}_1.fastq.gz", "out/{sample}/{sample}_2.fastq.gz"], sample=SAMPLES)

rule bam_to_fastq:
	input:
		config["input_path"]+"/{sample}.bam",
	output:
		fastq_1="out/{sample}/{sample}_1.fastq.gz",
                fastq_2="out/{sample}/{sample}_2.fastq.gz"
	conda: 'envs/samtools.yml'
	shell:
		"""
			samtools fastq -c 6 -1 {output.fastq_1} -2 {output.fastq_2} {input}
		"""
  
