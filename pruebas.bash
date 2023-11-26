#!/bin/bash
# Set the project root directory and paths
root="/home/kenneth/Documents/HLA/HLA"

#Iterate over the files
for i in "${root}"/$1*1.fastq.gz; do
    #Extract the files basenames
    if [[ $i =~ ($1[[:digit:]]+)_(.+)\.fastq\.gz ]]; then
        bs_number="${BASH_REMATCH[1]}"
        echo "${bs_number}"
    else
        echo "fail"
    fi
done
