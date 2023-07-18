import os
import glob
import pysam
from snakemake.utils import min_version


min_version("7.25.3")

# here we need wildcard contraint based on input file name pattern or it picks up bam file from sub-dirs
SAMPLES, = glob_wildcards(config["input_path"]+"/{sample}.Aligned.sortedByCoord.out.patched.md.bam")


FIRST_SAMPLE = str(SAMPLES[0])

def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in rule run_irap
    attemps = reiterations + 1
    Max number attemps = 6
    """
    mem_avail = [ 8, 16, 32, 64, 128, 300 ]  
    if attempt > len(mem_avail):
        print(f"Attemps {attempt} exceeds the maximum number of attemps: {len(mem_avail)}")
        print(f"modify value of --restart-times or adjust mem_avail resources accordingly")
        sys.exit(1)
    else:
        return mem_avail[attempt-1] * 1000


def detect_read_type(wildcards):
    if os.path.isfile(config["input_path"]+ f"/{wildcards['sample']}.Aligned.sortedByCoord.out.patched.md.bam.bai"):
        print(f"BAM index exists")
    else:
        print(f"generating index...")
        pysam.index(config["input_path"]+ f"/{wildcards['sample']}.Aligned.sortedByCoord.out.patched.md.bam")

    with pysam.AlignmentFile(config["input_path"]+ f"/{wildcards['sample']}.Aligned.sortedByCoord.out.patched.md.bam", "rb") as bam:
        for read in bam.fetch():
            if read.is_paired:
                return "pe"
    return "se"


def detect_read_type_first_sample(bam_file_name):
    if os.path.isfile(config["input_path"]+ f"/{bam_file_name}.Aligned.sortedByCoord.out.patched.md.bam.bai"):
        print(f"BAM index exists")
    else:
        print(f"generating index...")
        pysam.index(config["input_path"]+ f"/{bam_file_name}.Aligned.sortedByCoord.out.patched.md.bam")

    with pysam.AlignmentFile(config["input_path"]+ f"/{bam_file_name}.Aligned.sortedByCoord.out.patched.md.bam", "rb") as bam:
        for read in bam.fetch():
            if read.is_paired:
                return "pe"
    return "se"


rule all:
    input: expand(["out/{sample}/{sample}.fastq.val", "out/{sample}/{sample}_prepare_aggregation.done"], sample=SAMPLES), "out/workflow.done", f"out/stage0_{FIRST_SAMPLE}.txt"


rule check_bam:
    input:
        config["input_path"]+ "/{sample}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        bam_check = temp("out/{sample}/{sample}_1.bam_checked")
    log: "logs/{sample}/{sample}_check_bam.log"
    conda:
        "envs/samtools.yml"
    threads: 4
    params:
        bai = config["input_path"]+ "/{sample}.Aligned.sortedByCoord.out.patched.md.bam.bai"
    resources: 
        load=5
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        # ensure that the index file exists, otherwise create it
        if [ -e {params.bai} ]; then
            echo "BAM index file exists for {wildcards.sample}" 
        else
            echo "BAM index file does not exist for {wildcards.sample}. Generating it..." 
            samtools index -b --threads {threads} {input}
        fi

        samtools quickcheck {input}
        if [ $? -ne 0 ]; then
            echo "ERROR: {input} is not a valid BAM file"
        else
            echo "BAM file is valid"
            touch {output.bam_check}
        fi
        """

rule bam_to_fastq:
    """
    Produces interleaved Fastq from a sorted bam file
    """
    input:
        bam = config["input_path"]+"/{sample}.Aligned.sortedByCoord.out.patched.md.bam",
        check_bam = rules.check_bam.output.bam_check
    output:
        fastq = temp("out/{sample}/{sample}.fq"),
        sorted_bam= temp("out/{sample}/sorted_{sample}.bam")
    conda: 
        "envs/samtools.yml"
    threads: 8
    params: 
        read_type = detect_read_type
    log: "logs/{sample}/{sample}_bam_to_fastq.log"
    resources: 
        mem_mb=get_mem_mb,
        load=5
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        echo "read_type: {params.read_type}"

        samtools sort -n --threads {threads} -o {output.sorted_bam} {input.bam}

        # detect SE or PE
        if [[ {params.read_type} == "pe" ]]; then
            #If the input contains read-pairs which are to be interleaved or written to separate 
            #files in the same order, then the input should be first collated by name. Use samtools 
            #collate or samtools sort -n to ensure this.
            samtools fastq --threads {threads} -0 /dev/null -s /dev/null {output.sorted_bam} > {output.fastq}
        elif [[ {params.read_type} == "se" ]]; then
            samtools fastq --threads {threads} {output.sorted_bam} > {output.fastq}
        else
            echo "ERROR: read type not detected"
            exit 1
        fi
        """

	
checkpoint validating_fastq:
    """
    Here if else statement could be modified to run iRAP/ISL for PE and SE data.
    """
    input:
        fastq = rules.bam_to_fastq.output.fastq
    output:
        val_fastq = temp("out/{sample}/{sample}.fastq.val")
    log: "logs/{sample}/{sample}_validating_fastq.log"
    params:
        fastq = "out/{sample}/{sample}.fastq"
    conda:
        "envs/fastq_utils.yml"
    resources: 
        load=5
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
	
        fastq_info {input.fastq}  
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed fastq validation {input.fastq}"
            mv {input.fastq} {params.fastq}
        else
            echo "validation successful"
        fi

        touch {output.val_fastq}
        """

#def aggregate_fastq(wildcards):
#    checkpoint_output = checkpoints.validating_fastq.get(**wildcards).output[0]
#    print(checkoint_output)
#    return expand("out/{sample}/{fq}_fastqc.html",
#        sample=wildcards.sample,
#        fq=glob_wildcards(os.path.join(checkpoint_output, "{fq}.fastq.gz")).fq)

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
        irapMem=16000000000,
        irapDataOption="",
        filename= f"{FIRST_SAMPLE}",
        read_type = detect_read_type_first_sample(FIRST_SAMPLE)
    resources:
        mem_mb=16000,
        load=10
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

        pushd {params.root_dir}/scripts > /dev/null
        wget https://github.com/biopet/validatefastq/releases/download/v0.1.1/validatefastq-assembly-0.1.1.jar
        chmod +x validatefastq-assembly-0.1.1.jar
	    popd

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
 
        pushd $workingDir > /dev/null
	
        if [[ {params.read_type} == "se" ]]; then
            # fastq is SE
            cp {params.root_dir}/{input.fastq} $workingDir/${{localFastqPath}}.fastq
            java -jar {params.root_dir}/scripts/validatefastq-assembly-0.1.1.jar --fastq1 $workingDir/${{localFastqPath}}.fastq

            echo "Calling irap_single_lib...SE mode"
            cmd="irap_single_lib -0 -A -f -o irap_single_lib -1 $workingDir/${{localFastqPath}}.fastq -c {params.conf} -s {params.strand} -m {params.irapMem} -t {threads} -C {params.irapDataOption}"
            echo "stage0 will run now:"
            eval $cmd
            echo "stage0 finished"
        else
            # fastq is PE
            #split_fastq {input.fastq} $workingDir ${{localFastqPath}}
            reformat.sh ow=t int=t vpair=t vint=t in={params.root_dir}/{input.fastq} out1=$workingDir/${{localFastqPath}}_1.fastq out2=$workingDir/${{localFastqPath}}_2.fastq
            java -jar {params.root_dir}/scripts/validatefastq-assembly-0.1.1.jar --fastq1 $workingDir/${{localFastqPath}}_1.fastq --fastq2 $workingDir/${{localFastqPath}}_2.fastq

            echo "Calling irap_single_lib...PE mode"
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
    output: "out/{sample}/{sample}_irap_completed.done"
    conda: "envs/isl.yaml"
    log: "logs/{sample}/{sample}_irap.log"
    params:
        private_script=config["private_script"],
        conf=config["irap_config"],
        root_dir=config["atlas_gtex_root"],
        strand="both",
        irapDataOption="",
        filename="{sample}",
        first_sample=f"{FIRST_SAMPLE}",
        read_type = detect_read_type
    resources: 
        mem_mb=get_mem_mb,
        load=10
    threads: 16
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        irapMem=$(("{resources.mem_mb}000000"))
		
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
        fi

        pushd $workingDir > /dev/null
	

        if [[ {params.read_type} == "se" ]]; then
            # fastq is SE
            cp {params.root_dir}/{input.fastq} $workingDir/${{localFastqPath}}.fastq
            java -jar {params.root_dir}/scripts/validatefastq-assembly-0.1.1.jar --fastq1 $workingDir/${{localFastqPath}}.fastq

            echo "Calling irap_single_lib...SE mode"
            cmd="irap_single_lib -A -f -o irap_single_lib -1 $workingDir/${{localFastqPath}}.fastq -c {params.conf} -s {params.strand} -m $irapMem -t {threads} -C {params.irapDataOption}"
            echo "SE IRAP will run now:"
            eval $cmd
            echo "irap_single_lib SE finished for {wildcards.sample}"
        else
            # fastq is PE
            reformat.sh ow=t int=t vpair=t vint=t in={params.root_dir}/{input.fastq} out1=$workingDir/${{localFastqPath}}_1.fastq out2=$workingDir/${{localFastqPath}}_2.fastq
            java -jar {params.root_dir}/scripts/validatefastq-assembly-0.1.1.jar --fastq1 $workingDir/${{localFastqPath}}_1.fastq --fastq2 $workingDir/${{localFastqPath}}_2.fastq
            echo "Calling irap_single_lib..."
	    
            cmd="irap_single_lib -A -f -o irap_single_lib -1 ${{localFastqPath}}_1.fastq -2 ${{localFastqPath}}_2.fastq -c {params.conf} -s {params.strand} -m $irapMem -t {threads} -C {params.irapDataOption}"
            echo "PE IRAP will run now:"
            eval $cmd
            echo "irap_single_lib PE finished for {wildcards.sample}"
        fi

        popd

        # The files to collected by aggregation from irap are in $IRAP_SINGLE_LIB/out
        touch {output}
        """



rule prepare_aggregation:
    """
    Prepare aggregation of irap outputs, which are at
    $ISL_WORKING_DIR/irap_single_lib/{sample[0:5]}/{sample}
    and copy essential aggregation files into:
    $IRAP_SINGLE_LIB/out/{sample[0:5]}/{sample}
    """
    input: "out/{sample}/{sample}_irap_completed.done"
    output: "out/{sample}/{sample}_prepare_aggregation.done"
    log: "logs/{sample}/{sample}_prepare_agreggation.log"
    params:
        private_script=config["private_script"]
    resources: 
        load=2
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {params.private_script}/gtex_bulk_env.sh
        prefix_sample=$(echo {wildcards.sample} | cut -c1-6)

        # move irap outputs to $ISL_WORKING_DIR/irap_single_lib/out/sample[0:5]/sample

        destination_dir=$ISL_RESULTS_DIR/$prefix_sample/{wildcards.sample}
        mkdir -p $destination_dir

        cp $ISL_WORKING_DIR/irap_single_lib/$prefix_sample/{wildcards.sample}/{wildcards.sample}*kallisto*  $destination_dir/
        cp $ISL_WORKING_DIR/irap_single_lib/$prefix_sample/{wildcards.sample}/{wildcards.sample}*htseq2*  $destination_dir/
        cp $ISL_WORKING_DIR/irap_single_lib/$prefix_sample/{wildcards.sample}/{wildcards.sample}.versions.tsv  $destination_dir/irap.versions.tsv
        cp -r $ISL_WORKING_DIR/irap_single_lib/$prefix_sample/{wildcards.sample}/logs  $destination_dir/
        cp -r $ISL_WORKING_DIR/irap_single_lib/$prefix_sample/{wildcards.sample}/qc  $destination_dir/

        # file checks
        if ls "$destination_dir"/*kallisto* 1> /dev/null 2>&1; then
            echo "Files matching the pattern '*kallisto*' exist in the folder $destination_dir"
        else
            echo "No files matching the pattern '*kallisto*' exist in the folder $destination_dir"
            exit 1
        fi

        if ls "$destination_dir"/*htseq2* 1> /dev/null 2>&1; then
            echo "Files matching the pattern '*htseq2*' exist in the folder $destination_dir"
        else
            echo "No files matching the pattern '*htseq2*' exist in the folder $destination_dir"
            exit 1
        fi

        echo "tophat2 align summary:"
        cat $ISL_WORKING_DIR/processing_data/h/homo_sapiens/irap_qc/tophat2/*/{wildcards.sample}/{wildcards.sample}/align_summary.txt

        echo "removing other temp files in $ISL_RAW_DIR and $ISL_WORKING_DIR "
        rm -f $ISL_RAW_DIR/$prefix_sample/{wildcards.sample}/{wildcards.sample}*
        rm -rf $ISL_WORKING_DIR/irap_single_lib/$prefix_sample/{wildcards.sample}
        rm -rf $ISL_WORKING_DIR/processing_data/h/homo_sapiens/irap_qc/tophat2/*/{wildcards.sample}
        rm -rf $ISL_WORKING_DIR/processing_data/h/homo_sapiens/irap_qc/none/kallisto/*/{wildcards.sample}
        rm -rf $ISL_WORKING_DIR/processing_data/h/homo_sapiens/irap_qc/tophat2/htseq2/*/{wildcards.sample}

        echo "GTEX analysis completed for {wildcards.sample}. Ready for aggregation"

        touch {output}
        """


#rule aggregate_libraries:
#    """
#    Final rule to agregate all library outputs, and store them in a single folder
#    $IRAP_SINGLE_LIB/studies/E-GTEX-8/homo_sapiens/
#    """



rule final_workflow_check:
    input: expand(["out/{sample}/{sample}_prepare_aggregation.done"], sample=SAMPLES)
    output: "out/workflow.done"
    resources: 
        load=1
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        touch {output}
        """
