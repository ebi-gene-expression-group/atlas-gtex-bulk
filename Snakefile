import os
import glob

SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.bam")

# probably need to specify here (or in config) the species or the location of the genome/annotations

rule all:
    input: expand(["out/{sample}/{sample}_1_fastqc.html","out/{sample}/{sample}_2_fastqc.html","out/{sample}/{sample}_1.fastq.gz", "out/{sample}/{sample}_2.fastq.gz"], sample=SAMPLES)


rule bam_to_fastq:
	input:
		bam = config["input_path"]+"/{sample}.bam",
	output:
		fastq_1 = "out/{sample}/{sample}_1.fq.gz",
                fastq_2 = "out/{sample}/{sample}_2.fq.gz"
	conda: 
		"envs/samtools.yml"
	shell:
		"""
			samtools fastq -c 6 -1 {output.fastq_1} -2 {output.fastq_2} {input.bam}
		"""
 
rule validating_fastq:
        input:
                fastq_1 = rules.bam_to_fastq.output.fastq_1,
                fastq_2 = rules.bam_to_fastq.output.fastq_2,
        output:
                val_fastq = "out/{sample}/{sample}_1.fastq.val",
		fastq_1 = "out/{sample}/{sample}_1.fastq.gz",
                fastq_2 = "out/{sample}/{sample}_2.fastq.gz",
        conda:
                "envs/fastq_utils.yml"
	shell:
                """
                	fastq_info {input.fastq_1} {input.fastq_2} 
			if [ $? -ne 0 ]; then
                        	rm -rf {input.fastq_1}
                        	rm -rf {input.fastq_2}
                        	echo "ERROR: Failed fastq validation {input.fastq_1} and {input.fastq_2}"
                    	else
				echo "validation successful"
                    		mv {input.fastq_1} {output.fastq_1}
                    		mv {input.fastq_2} {output.fastq_2}
			fi

			touch {output.val_fastq}
		"""


rule qc:
	input:
		fq1 = rules.validating_fastq.output.fastq_1,
		fq2 = rules.validating_fastq.output.fastq_2,
		check = rules.validating_fastq.output.val_fastq
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


