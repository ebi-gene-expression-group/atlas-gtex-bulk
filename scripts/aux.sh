#!/usr/bin/env bash

get_local_relative_library_path()
{
	local library=$1;
      	local subDir=${2:-''};
     	local remoteLibraryPath=`get_library_path $library $subDir`;
      	local localLibraryPath="raw_data/$subDir"`echo $remoteLibraryPath`;
     	echo $localLibraryPath
}

get_library_path()
{
    local library=$1;
    local rootDir=${2:-''};
    local forceShortForm=${3:-''};
    local subDir=${library:0:6};
    #local prefix=;
    #if ! [[ $subDir =~ "ENC" ]] && [[ -z "$forceShortForm" ]]; then
    #    local num=${library:3};
    #    if [ $num -gt 1000000 ]; then
    #        digits=$(echo ${library:9} | sed 's/^0*//');
    #        prefix="$(printf %03d $digits)/";
    #    fi;
    #fi;
    echo "${rootDir}${subDir}/${library}/${library}"
}

split_fastq()
{
    # some checks could be implemented here
    local fastq=$1;
    local workdir=$2;
    local pathToFile=$3;

    cat $fastq | grep -E '^@[^\s]+ 1[^\n]+$|^@[^\/\s]+\/1+$' -A 3 --no-group-separator > $workdir/${pathToFile}_1.fastq
    cat $fastq | grep -E '^@[^\s]+ 2[^\n]+$|^@[^\/\s]+\/2+$' -A 3 --no-group-separator > $workdir/${pathToFile}_2.fastq
    # altertive regex:
    # ^@[^\s]+ 1[^\n]+$|^@[^\/\n]+\/1+$ --> this is for read 1
    # ^@[^\s]+ 2[^\n]+$|^@[^\/\n]+\/2+$ --> this is for read 2    
    echo "split fastq completed"
}
