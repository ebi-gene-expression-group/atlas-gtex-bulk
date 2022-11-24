import os

SAMPLES = ["mt.sorted"]

rule all:
    input: expand(["{sample}_1.fastq", "{sample}_2.fastq"], sample=SAMPLES)

rule bam_to_fastq:
	input:
		"test-data/{sample}.bam",
	output:
		fastq_1="{sample}_1.fastq",
                fastq_2="{sample}_2.fastq"
	conda: 'envs/samtools.yml'
	shell:
		"""
			samtools fastq -F 2816 -c 6 -1 {output.fastq_1} -2 {output.fastq_2} {input}
		"""
  
