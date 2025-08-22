#!/bin/bash

mkdir -p crispr_site

for fasta_file in $(ls $1); do
   
    exomename=$(echo "$fasta_file" | awk -F'_' '{print $1}')
    grep -E -B1   '.{20,}.GG' $1/"$fasta_file" | grep -v "^--" > "crispr_site/${exomename}_precrispr.fasta"

done       


###cml ./identifyCrisprSite.sh topmotifs
