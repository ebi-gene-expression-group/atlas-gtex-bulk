#!/usr/bin/env bash

source $IRAP_SINGLE_LIB/gtex_bulk_env.sh

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${scriptDir}/../isl/lib/irap.sh

# Input if paired-end or single-end
sepe="PAIRED"

# Add logic here to get localFastqPath
localFastqPath=""

# Add logic here to get conf
conf=""

# Add logic here to get strand
# Get from AE2_BASE_DIR sdrf, if "not applicable" then replace it with "both"
# cat /nfs/production/irene/ma/ae2_production/data/EXPERIMENT/GTEX/E-GTEX-8/E-GTEX-8.sdrf.txt | awk -F "\t" '{ print $29 }'
strand=""

# Add logic here to get mem
lsfMem=4096 ### This is the minimum, set when STUDIES.LSF_MEM is <null>
### if failing due to insufficient mem, increase this (insert logic)
irapMem=$(($lsfMem*1000000)) 

irapDataOption="-i data_dir=$workingDir/data"

# Execute below if a sepe logic is ready
# if [ "$sepe" == "PAIRED" ]; then
#     cmd="irap_single_lib$tidy -A -f -o irap_single_lib -1 ${localFastqPath}_1.fastq.gz -2 ${localFastqPath}_2.fastq.gz -c $conf -s $strand -m $irapMem -t 5 -C $irapDataOption"
# else
#     cmd="irap_single_lib$tidy -A -f -o irap_single_lib -1 ${localFastqPath}.fastq.gz -c $conf -s $strand -m $irapMem -t 5 -C $irapDataOption"
# fi 



# Below are tests

# ========== 
# Main test command
# ==========

localFastqPath="/homes/irisyu/gtex_bulk/tests/test_data/GSM461177"
conf=$ISL_CONFIG_DIR/homo_sapiens.conf
strand="both"
lsfMem=4096
irapMem=$(($lsfMem*1000000))
workingDir=$ISL_WORKING_DIR
irapDataOption="-i data_dir=$workingDir/data"

# irap_single_lib -A -f -o irap_single_lib -1 ${localFastqPath}_1.fastqsanger.gz -2 ${localFastqPath}_2.fastqsanger.gz -c $conf -s $strand -m $irapMem -t 5 -C $irapDataOption


# ==========
# irap command to trace stage0 incomplete error
# ==========
cmd=irap
irap_options="-i data_dir=$workingDir/data isl_mode=y"
DATA_DIR=`conf_get_data_dir $conf $irap_options`
# for conf_get_data_dir, see https://github.com/nunofonseca/irap/blob/4e04521f3cd00718cb3ea998d7ef81025c484117/aux/sh/irap_shared.sh#L30
mem=$irapMem   ### same as above
threads=5

irap conf=$conf data_dir=$DATA_DIR $irap_options max_mem=$mem max_threads=$threads stage0 -n -q
# =====