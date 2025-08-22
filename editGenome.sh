#!/bin/bash

# Make Directory to save results 
mkdir -p edit_Genome

# for loop to insert A before NGG 
for fasta_file in $(ls $1); do
        
        exomename=$(echo "${fasta_file}" | awk -F'_' '{print $1}')
        
        sed '/^>/!s/\(.GG\)/A\1/g'  $1/$fasta_file | grep -v "^--" > "edit_Genome/${exomename}_postcrispr.fasta"

done 
 
##cml ./bash editGenome.sh  crispr_site
