#!/bin/bash

# $1 and $2 argumets are corresponding to exomesCohort (contains fasta files)  folder and motif_list.txt (contains list of motifs)

# create topmotifs directory to save *_topmotifs.fasta files

mkdir -p  topmotifs


#for loop over each exome inside the exomesCohort  folder  

for fasta in $(ls "${1}");do
  
    #make temporary text file

    touch tm_file.txt
    
    #nested for loop to calculate the  occurrance of all  motifs in each exome    

    for pattern in $(cat "${2}");do
    
        count=$(awk -v motif="${pattern}" 'BEGIN { count=0 } { count += gsub(motif, "", $0) } END { print count }' "${1}/${fasta}")

        echo "${pattern}  ${count}" >> tm_file.txt

    done
 
    #identify the 3 highest occurring motifs in each exome and save it in temporary txt file

    sort -k2,2nr  tm_file.txt | head -n 3     | awk '{print $1}' > "top_occurance.txt"

    # creat the out out file name for each exome    
    exomename=$(echo "$fasta" | awk -F'.' '{print $1}')

    # dynamically build grep pattern for top motifs
    top_motifs=$(paste -sd '|' top_occurance.txt)
    
    # Output the genes with 3 top motifs and corresponding sequences to a new file
    grep -B1 -E "$top_motifs" "$1/${fasta}" | grep  -v "^--" | awk '!seen[$0]++' >> "topmotifs/${exomename}_topmotifs.fasta"
     
    # remove temporary text files
    rm tm_file.txt 
    rm top_occurance.txt

done


   

