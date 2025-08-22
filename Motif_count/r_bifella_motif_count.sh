#!/bin/bash
#week_2 assignment
#Author: Azadeh Nikouee

#create motifs directory
mkdir -p motifs

#create motif_count.txt file
touch motif_count.txt

#creat list of motifs 
motifs=$(cat interesting_motifs.txt)

#create for loop to count genes with target motif and count them
for pattern in $motifs ;
    do
    #count motif in fasta file
    motif_count=$(grep -c "$pattern" r_bifella.fasta) 
    echo "$pattern $motif_count" >> motif_count.txt

    #extract genes and sequence and save each in motifs directory
    grep -B 1 $pattern r_bifella.fasta > ./motifs/${pattern}.fasta

done

echo "file creation completed"



