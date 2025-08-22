#!/bin/bash

## $1=clinical_data.txt   
## $2=exomes (direcory of FASTA files)

#replace  space between names with two parts by _ and organize columns and then filter file by awk for column 3 and column 5 and extract code_name(column 6)
 
i() { 
 tr ' ' '_' < "$1" | tr '\t' ' ' | column -t |
  awk '$3 >= 20 && $3 <= 30 && $5 == "Sequenced" {print $6} ' 
 }


# make directory

mkdir -p  exomesCohort


# Loop through filtered sample names and copy matching FASTA files in exomes
for sample in $(i "$1"); do
    fasta_file="$2/${sample}.fasta"  # Path to FASTA file

    cp $fasta_file exomesCohort

done   
