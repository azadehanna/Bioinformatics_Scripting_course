import requests  # make API requests
import json  #for JSON responses
import re  #search for start codons in DNA sequences
import os  # for file and directory 
from Bio.Seq import Seq  # for translation
import sys

#  get Ensembl Gene ID (MC1R)

MYGENE_URL = "https://mygene.info/v3/query?q=MC1R&species=human&fields=ensembl.gene"

response = requests.get(MYGENE_URL)
data = response.json()
print(data)
ensembl_id = data['hits'][0]['ensembl']['gene']
print(f"Ensembl ID for MC1R: {ensembl_id}")

#  get the Nucleotide Sequence of MC1R

ENSEMBL_URL = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=application/json"
response = requests.get(ENSEMBL_URL)
sequence_data = response.json()
nucleotide_seq = sequence_data['seq']
#print (nucleotide_seq)

# find the longest open reading frame 
def find_longest_orf(dna_seq):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    longest_orf = ""

    for frame in range(3):
        seq = dna_seq[frame:]
        start_positions = [m.start() for m in re.finditer(start_codon, seq)]
        
        for start in start_positions:
            for i in range(start, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if codon in stop_codons:
                    orf = seq[start:i+3]
                    if len(orf) > len(longest_orf):
                        longest_orf = orf
                    break
    
    return longest_orf

longest_orf = find_longest_orf(nucleotide_seq)

protein_seq = Seq(longest_orf).translate()

# write the Sequences to FASTA File

fasta_path = os.path.join("./", "mc1r_sequence.fasta")

with open(fasta_path, "w") as fasta_file:
    fasta_file.write(f">MC1R_Nucleotide_Sequence\n{nucleotide_seq}\n")
    fasta_file.write(f">MC1R_Longest_ORF_Protein\n{protein_seq}\n")

#print(f"FASTA file created: {fasta_path}")


server = "https://rest.ensembl.org"
ext = f"/homology/id/human/{ensembl_id}?"

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

decoded = r.json()

species_set = set()  # Use a set to remove duplicates automatically

# Extract species from homology data
for entry in decoded.get("data", []):  # Ensure "data" exists
    for homology in entry.get("homologies", []):  # Ensure "homologies" exists
        species = homology.get("target", {}).get("species")
        if species:
            species_set.add(species)  # Add to set to avoid duplicates

# Convert set to sorted list
unique_species_list = sorted(species_set)

# Save to a text file
with open("mc1r_homology_list.txt", "w") as file:
    for species in unique_species_list:
        file.write(species + "\n")

print("Extracted species names saved in mc1r_homology_list.txt")

 
