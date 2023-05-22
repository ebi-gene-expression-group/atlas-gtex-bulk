import os
import glob
from snakemake.utils import min_version


min_version("7.25.3")

# here we need wildcard contraint based on input file name pattern or it picks up bam file from sub-dirs
SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.Aligned.sortedByCoord.out.patched.md.bam")

# probably need to specify here (or in config) the species or the location of the genome/annotations
FIRST_SAMPLE = str(SAMPLES[0])

def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules 
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [ 2, 2, 4, 8, 16, 64, 128, 300 ]  
    if attempt > len(mem_avail):
        print(f"Attemps {attempt} exceeds the maximum number of attemps: {len(mem_avail)}")
        print(f"modify value of --restart-times or adjust mem_avail resources accordingly")
        sys.exit(1)
    else:
        return mem_avail[attempt-1] * 1000




rule all:
    input: expand(["out/{sample}/{sample}.fastq.val", "out/{sample}.txt"], sample=SAMPLES), "done.txt", f"out/stage0_{FIRST_SAMPLE}.txt"


rule check_bam:
    input:
        config["input_path"]+ "/{sample}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        bam_check = "out/{sample}/{sample}_1.bam_checked"
    conda:
        "envs/samtools.yml"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set

        samtools quickcheck {input}
        if [ $? -ne 0 ]; then
            echo "ERROR: {input} is not a valid BAM file" >&2
        else
            touch {output.bam_check}
        fi
        """

rule bam_to_fastq:
    input:
        bam = config["input_path"]+"/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        check_bam = rules.check_bam.output.bam_check
    output:
        fastq = "out/{sample}/{sample}.fq",
        sorted_bam="out/{sample}/sorted_{sample}.bam"
    conda: 
        "envs/samtools.yml"
    threads: 8
    resources: mem_mb=8000
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        samtools sort -n --threads {threads} -o {output.sorted_bam} {input.bam}

        #If the input contains read-pairs which are to be interleaved or written to separate 
        #files in the same order, then the input should be first collated by name. Use samtools 
        #collate or samtools sort -n to ensure this.
        samtools fastq --threads {threads} -0 /dev/null {output.sorted_bam} > {output.fastq}
        """

	
checkpoint validating_fastq:
    """
    Here if else statement could be modified to run iRAP/ISL for PE and SE data.
    """
    input:
        fastq = rules.bam_to_fastq.output.fastq
    output:
        val_fastq = "out/{sample}/{sample}.fastq.val"
    params:
        fastq = "out/{sample}/{sample}.fastq"
    conda:
        "envs/fastq_utils.yml"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
	
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

rule fastqc:
    """
    Not currently used
    """
    input:
        check = rules.validating_fastq.output.val_fastq
    output:
        fqc = "out/{sample}/{fq}_1_fastqc.html"
    params:
        fq = "out/{sample}/{fq}.fastq",
        fqc_dir = "out/{sample}"
    conda:
        "envs/fastqc.yml"
    threads: 4
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        fastqc -c {threads} {params.fq} --outdir={params.fqc_dir}
        """

rule run_irap_stage0:
    """
    This ensures Irap stage0 is run only once, for the first sample
    """
    input:
        fastq= f"out/{FIRST_SAMPLE}/{FIRST_SAMPLE}.fq",
        check = f"out/{FIRST_SAMPLE}/{FIRST_SAMPLE}.fastq.val"
    output: f"out/stage0_{FIRST_SAMPLE}.txt"
    conda: "envs/isl.yaml"
    log: f"logs/irap_stage0_{FIRST_SAMPLE}.log"
    params:
        private_script=config["private_script"],
        conf=config["irap_config"],
        root_dir=config["atlas_gtex_root"],
        strand="both",
        irapMem=4096000000,
        irapDataOption="",
        filename= f"{FIRST_SAMPLE}"
    resources: mem_mb=10000
    threads: 16
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
		
        source {params.private_script}/gtex_bulk_env.sh
        source {params.private_script}/gtex_bulk_init.sh
        #source {params.root_dir}/isl/lib/functions.sh
        source {params.root_dir}/isl/lib/generic_routines.sh
        source {params.root_dir}/isl/lib/process_routines.sh
        source {params.root_dir}/isl/lib/irap.sh
	
        cp {params.private_script}/gtex_bulk_env.sh $IRAP_SINGLE_LIB
	
        cat {input.check}

        library={params.filename}
        echo "library: $library"
        workingDir=$ISL_WORKING_DIR

        source {params.root_dir}/scripts/aux.sh

        which get_local_relative_library_path
        which get_library_path


        localFastqPath=$(get_local_relative_library_path $library )
	
        echo "workingDir: $workingDir"
        echo "localFastqPath: $localFastqPath"
        
        mkdir -p $(dirname $workingDir/$localFastqPath)

        #split_fastq {input.fastq} $workingDir ${{localFastqPath}}
        reformat.sh ow=t int=t vpair=t vint=t in={input.fastq} out1=$workingDir/${{localFastqPath}}_1.fastq out2=$workingDir/${{localFastqPath}}_2.fastq

        pushd $workingDir > /dev/null
	
        sepe=$( fastq_info $workingDir/${{localFastqPath}}_1.fastq $workingDir/${{localFastqPath}}_2.fastq )

        if [ $? -ne 0 ]; then
            # fastq is SE
            # iRAP SE command here
            echo "SE "
        else
            # fastq is PE
            echo "Calling irap_single_lib..."
            cmd="irap_single_lib -0 -A -f -o irap_single_lib -1 ${{localFastqPath}}_1.fastq -2 ${{localFastqPath}}_2.fastq -c {params.conf} -s {params.strand} -m {params.irapMem} -t {threads} -C {params.irapDataOption}"
            echo "stage0 will run now:"
            eval $cmd
            echo "stage0 finished"
        fi

        popd

        touch {output}
        """


rule run_irap:
    input:
        fastq=rules.bam_to_fastq.output.fastq,
        check = rules.validating_fastq.output.val_fastq,
        stage0_completed=rules.run_irap_stage0.output
    output: "out/{sample}.txt"
    conda: "envs/isl.yaml"
    log: "logs/{sample}_irap.log"
    params:
        private_script=config["private_script"],
        conf=config["irap_config"],
        root_dir=config["atlas_gtex_root"],
        strand="both",
        irapMem=4096000000,
        irapDataOption="",
        filename="{sample}",
        first_sample=FIRST_SAMPLE
    resources: mem_mb=10000
    threads: 16
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
		
        source {params.private_script}/gtex_bulk_env.sh
        #source {params.private_script}/gtex_bulk_init.sh
        #source {params.root_dir}/isl/lib/functions.sh
        source {params.root_dir}/isl/lib/generic_routines.sh
        source {params.root_dir}/isl/lib/process_routines.sh
        source {params.root_dir}/isl/lib/irap.sh
	
        ##cp {params.private_script}/gtex_bulk_env.sh $IRAP_SINGLE_LIB
	
        cat {input.check}

        library={params.filename}
        echo "library: $library"
        workingDir=$ISL_WORKING_DIR

        source {params.root_dir}/scripts/aux.sh

        which get_local_relative_library_path
        which get_library_path


        localFastqPath=$(get_local_relative_library_path $library )
	
        echo "workingDir: $workingDir"
        echo "localFastqPath: $localFastqPath"

        echo "sample: {wildcards.sample}"
        echo "first sample: {params.first_sample}"

        if [[ "{wildcards.sample}" != "{params.first_sample}" ]]; then
            mkdir -p $(dirname $workingDir/$localFastqPath)
            #split_fastq {input.fastq} $workingDir ${{localFastqPath}}
            reformat.sh ow=t int=t vpair=t vint=t in={input.fastq} out1=$workingDir/${{localFastqPath}}_1.fastq out2=$workingDir/${{localFastqPath}}_2.fastq
        fi

        pushd $workingDir > /dev/null
	
        sepe=$( fastq_info $workingDir/${{localFastqPath}}_1.fastq $workingDir/${{localFastqPath}}_2.fastq )

        if [ $? -ne 0 ]; then
            # fastq is SE
            # iRAP SE command here
            echo "SE "
        else
            # fastq is PE
            echo "Calling irap_single_lib..."
	    
            cmd="irap_single_lib -A -f -o irap_single_lib -1 ${{localFastqPath}}_1.fastq -2 ${{localFastqPath}}_2.fastq -c {params.conf} -s {params.strand} -m {params.irapMem} -t {threads} -C {params.irapDataOption}"
            echo "PE IRAP will run now:"
            eval $cmd
            echo "irap_single_lib finished for {wildcards.sample}"
        fi

        popd

        # The files to collected by aggregation from irap are in $IRAP_SINGLE_LIB/out
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

rule merge:
    input: expand(["out/{sample}.txt"], sample=SAMPLES)
    output: "done.txt"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        touch {output}
        """
