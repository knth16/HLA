#!/bin/bash

################################################################################
# Script: hlahd_analysis_script.sh
# Author: Kenneth Valerio Aguilar
# Date: 11/18/2023
# Version: 1.0
# Purpose: Perform HLA genotyping analysis using hlahd.sh on multiple individuals
# Input Requirements: Requires fastq files with specific naming convention
# Usage: ./hlahd_analysis_script.sh
# Output: Log files, result files, and summary table in the results directory
################################################################################

# Set the project root directory
root="/molbio/projects/hla_genotyping/hlahd"

#Iterate over the files
for i in $root/fastq/BS*1.fastq.gz; do
    #Extract the files basenames
    name=basename ${i}1.fastq.gz 
    #Uncrompres the files into fastqtemp directory
    gunzip -c name1.fastq.gz  > fastqtmp


if [[ $filepath =~ BS([[:digit:]]+)_(.+)\.fastq\.gz ]]; then
    bs_number="${BASH_REMATCH[1]}"
    before_underscore="${BASH_REMATCH[2]%%_*}"