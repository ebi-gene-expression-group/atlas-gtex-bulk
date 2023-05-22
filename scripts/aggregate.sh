#!/usr/bin/env bash

# This script collates irap_single_lib results of individual libraries into one set of outputs for a given study.
# This is intended to run as a non-fg_atlas user, without access to the ISL database.
# The direct application of this script is to aggregate single-library results of a controlled-access study, such as GTEx.

source $IRAP_SINGLE_LIB/gtex_bulk_env.sh

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${scriptDir}/../isl/lib/functions.sh

studyId=$1
organism="homo_sapiens"
lsfMem=${2:-4096}

check_variables 'studyId'

# Have some sort of mechanism to log memory usage during aggregation,
# so that when aggregation fails, higher mem can be requested in the next run.
# In ISL, the memory record is in the db.  Here, we can use some a simple log.

# Remove any existing aggregation outputs

echo "Cleaning up for ${organism} in ${studyId}"
remove_process_working_dir "$studyId" aggregate "$organism"
remove_process_results "$studyId" aggregate "$organism"
echo "Done study cleanup"

# Generate a string listing directories with results to be aggregated. 
# Also copy irap.versions.tsv files to $workingDir - all of them should be
# the same for each library of the same study. If not, copy each different
# version of this file into irap.versions.tsv.<version> - for inspection
# later - and log there's more than one version of the file

workingDir=$(process_working_dir "$studyId" aggregate "$organism" 'yes')   ### correct?
versions="${workingDir}/irap.versions.tsv"
versionString=

libraryPathsForStudy=

# Try to get libraries from SDRF files first (preferred)

libs=$(get_libraries_for_ae2_experiment $studyId)

if [ $? -ne 0 ]; then

    # Getting libs from the config SDRF failed, so fall back to using the 
    # config XML file to get runs. This happens when we get crufty SDRF files,
    # which is sometimes be unavoidable for controlled-access experiments,
    # e.g., the expt is a mix of public (ENA) and controlled-access libraries.
    warn "Failed to derive runs from SDRF file for $studyId, fetching all libraries from the config XML file instead."
    libs=$(grep -P "<assay[^>]{0,}>[^<]+</assay>" $ATLAS_PROD/analysis/*/rna-seq/experiments/$studyId/$studyId-configuration.xml | awk -F '[<>]' '{ print $3 }')
    
    if [ $? -ne 0 ]; then
        warn "[ERROR] Can't get libraries to aggregate for ${studyId}-${organism}"
        return 1
    fi
fi

while read -r library; do
    resultsDir=$(process_results_dir $library irap_single_lib $organism no no)
    libraryPathsForStudy="$libraryPathsForStudy $resultsDir"

    # Check versions. Take out the parameters, which aren't helpful for
    # users anyway, specifying config file paths etc. Deal with either
    # versions.tsv or irap.vesions.tsv, whichever is available
 
    versionsFile="${resultsDir}/versions.tsv"
    if [ ! -s "$versionsFile" ]; then
        versionsFile="${resultsDir}/irap.versions.tsv"
    fi

    libVersions=$(cat "$versionsFile" | grep -v "iRAP params")
    
    if [ -z "$versionString" ]; then
        versionString=$libVersions
        echo -e "$versionString" >> $versions
    elif [ "$libVersions" != "$versionString" ]; then
        echo -e "$versionString" >> ${versions}.${library}
    fi

done <<< "$libs"

# Write library paths to file to pass to iRAP for aggregation

aux=${ISL_TEMP_DIR}/aux.${studyId}.${organism}.$$
echo $libraryPathsForStudy > $aux
            
# Submit aggregation to cluster

aggrCmd="irap_single_lib2report_atlas -B -j 4 folders_file=$aux out=$ISL_WORKING_DIR/studies/$studyId/$organism name=$studyId"
echo "${aggrCmd}"   ###
isl_lsf_submit 'aggregate' ${studyId}.${organism} "/irap_gtex" $lsfMem 1 $ISL_WORKING_DIR "$aggrCmd"

if [ $? -ne 0 ]; then
    echo  "[ERROR] Failed to submit to lsf job: '$aggrCmd'" 
    exit 1
else
    echo "Submitted aggregation for ${studyId} ${organism}"

    # Monitor aggregation job
    monitorLogPrefix=$(get_log_prefix processing 'aggregate' ${studyId}.${organism})
    monitorLog="${monitorLogPrefix}.log"
    isl_lsf_monitor 'aggregate' "${studyId}.${organism}" $monitorLog

    # Once done, move to the fg_atlas ISL studies dir.
    if [ $? -ne 0 ]; then
        echo "LSF processing completed with errors"
        exit 1
    else
        echo "LSF processing completed, no errors"

        # Move results to results directory
        resultsDir=$(process_results_dir $studyId aggregate $organism 'no')
        rm -rf $resultsDir && mkdir -p $(dirname $resultsDir) && mv $workingDir $resultsDir            
    fi
fi

# Once move is done, generate a file that indicates that the process is done.
# Below path can still be changed to a more logical location
touch $IRAP_SINGLE_LIB/$studyId.complete