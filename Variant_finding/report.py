
#!/usr/bin/env python3
## report  for result of pipeline 
import pandas as pd
import os
import csv
import ast

input_1 = "./harrington_clinical_data.txt"
input_2 = "./variant_results"
clinical_data = pd.read_table(input_1)

def read_text_files(input_2):
    """Reads all text files in a directory and returns their content as a list."""
    file_contents = []
    for filename in os.listdir(input_2):
        if filename.endswith(".txt"):
            filepath = os.path.join(input_2, filename)
            try:
                with open(filepath, 'r') as file:
                    for line in file:
                        list_1 = line.strip().split()
                    content = file.read()
                    file_contents.append(content)
            except Exception as e:
                 print(f"Error reading {filename}: {e}")
    return file_contents

# Initialize a list to store variant results for each sample.
variant_results = []

# Loop through each text file in the variant_results directory.
for filename in os.listdir(input_2):
    if filename.endswith(".txt"):
        sample_name = filename.split('_')[0]  # extract sample name from filename
        filepath = os.path.join(input_2, filename)
        mutation_found = False
        
        # Open and parse the file as a tab-delimited file.
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Check if a mutation was detected.
                if row["VariantDetected"].strip() == "True":
                    mutation_found = True
                    position = row["Position"]
                    coverage = row["Coverage"]
                    mutation_str = row["Mutation"].strip()
                    
                    # If mutation is in the format "G>C", extract the mutated base.
                    if ">" in mutation_str:
                        _, mutated_base = mutation_str.split(">")
                    else:
                        mutated_base = mutation_str
                        
                    # Convert the BaseFrequencies string into a dictionary.
                    base_freqs = ast.literal_eval(row["BaseFrequencies"])
                    frequency = base_freqs.get(mutated_base, "N/A")
                    
                    variant_results.append({
                        "Name": sample_name,
                        "position": position,
                        "coverage": coverage,
                        "mutation": mutated_base,
                        "frequency": frequency
                    })
                    break  # stop after finding the first mutation row
        # If no mutation was detected in this file, add a placeholder result.
        if not mutation_found:
            variant_results.append({
                "Name": sample_name,
                "position": None,
                "coverage": None,
                "mutation": None,
                "frequency": None
            })

# Convert the variant results list into a DataFrame.
variant_df = pd.DataFrame(variant_results)

# Merge the clinical data with the variant results on the Name column.
merged_data = pd.merge(clinical_data, variant_df, on="Name", how="left")

# Print out the report for each sample.
with open("report.txt", "w") as report_file:
    for idx, row in merged_data.iterrows():
        if row["mutation"]:
            report_file.write(
                f"Sample {row['Name']} had a mold color {row['Color']}, {row['coverage']} reads, and had {row['frequency']}% of the reads at position {row['position']} with the mutation {row['mutation']}.\n"
            )
        else:
            report_file.write(
                f"Sample {row['Name']} had a mold color {row['Color']} with no mutation detected.\n"
            )