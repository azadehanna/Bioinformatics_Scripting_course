import os
import csv

clinical_data_file = "clinical_data.txt"      #Input file       
report_filename = "exome_report.txt"          #output file
filtered_organisms = []

# read clinical data and filter codename
with open(clinical_data_file, "r") as f:
    lines = f.readlines()
    headers = lines[0].strip().split("\t")
    for line in lines[1:]:
            row = line.strip().split("\t")
            diameter = int(row[2])
            ## define cut off and condition for selected coloumn
            if 20 <= diameter <= 30 and row[4].lower() == "sequenced":
                filtered_organisms.append({
                    "discoverer": row[0],
                    "diameter": diameter,
                    "codename": row[5].strip().lower(),
                    "environment": row[3]
                })


### define function to read fasta file and select the first gene and it's sequence
def read_first_fasta(fasta_file):
        with open(fasta_file, "r") as file:
            header = None
            sequence = []
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if header:  # if we already have a header, stop at the next one
                        break
                    header = line  # capture the first header
                else:
                    sequence.append(line)  # capture sequence lines
            return header, ''.join(sequence) 



# Write the final report
with open(report_filename, "w") as report:
    for org in filtered_organisms:
        # the path of fasta files
        fasta_file_path = "./edit_Genome/" + org['codename'] + "_postcrispr.fasta"
        
        # read the first gene and sequence
        fasta_header, fasta_seq = read_first_fasta(fasta_file_path)

        # write the information of each organism from clinical data
        report.write("Organism " + org['codename'].capitalize() + ", discovered by " + org['discoverer'] + ", has a diameter of " + str(org['diameter']) + "mm, and is from the " + org['environment'] + " environment.\n")
        report.write("The list of genes can be found in: " + fasta_file_path + "\n")

        # write the first gene and it's sequence
        if fasta_header and fasta_seq:
            report.write("The first sequence of " + org['codename'] + " is:\n\n" + fasta_header + "\n" + fasta_seq + "\n\n")




