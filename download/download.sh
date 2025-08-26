#!/usr/bin/env bash

manifest_name=$1

if [ -z "$manifest_name" ]; then
    echo "Please supply the name of a manifest file you generated using instructions at https://anvilproject.org/learn/reference/gtex-v8-free-egress-instructions" 1>&2
    exit 1
fi

manifest="manifests/${manifest_name}.json"

if [ ! -e "$manifest" ]; then
    echo "No manifest at $manifest" 1>&2
    exit 1
else
    mkdir -p data/$manifest_name
fi

command="echo y | ./bin/gen3-client download-multiple --numparallel 10 --skip-completed --profile=gtex --manifest=$(pwd)/$manifest --download-path=$(pwd)/data/$manifest_name --protocol=s3"
echo $command
bsub -u jmanning -o $(pwd)/download.out.txt -e $(pwd)/download.err.txt "$command"
