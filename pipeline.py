#!/usr/bin/env python3
import os
import gzip
import subprocess
import pysam

##### Start  of Demultiplexing-section 1A and 1B

############### Start of ParseFastQ Class ###############
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self, filePath, headerSymbols=['@', '+']):
        ## check if the file is zipped by examining extension
        if filePath.endswith('.gz'):
            #open the zipped file in text mode("rt")
            self._file = gzip.open(filePath, 'rt')
        else:
            #open text file in readmode ("r")
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        #store the header symbolees(default:"@" for sequence header ,"+" for quality header)
        self._hdSyms = headerSymbols

    def __iter__(self):
        #make this object an interator by returning itself
        return self

    def __next__(self):
        # create an empty list to store the four parts of FastQ 
        elemList = []
        #loop four times to read each of four lines in a FastaQ
        for i in range(4):
            #read oneline from the file
            line = self._file.readline()
            #Increment the line counter regardless of whether a line was read
            self._currentLineNumber += 1
            if line:
                #If a line is read, remove the trailing newline and add it to elmList
                elemList.append(line.strip('\n'))
            else:
                #If no line is read(end-of-file), appened none to keep list length
                elemList.append(None)
        #If all four enteries in elemList are None, we've reached the end of the file
        if elemList.count(None) == 4:
            raise StopIteration
        # validate that no part of the record is missing
        assert all(elemList), "** ERROR: Incomplete FastQ record near line %s" % self._currentLineNumber
        #validate the first line(sequence header)starts with the expected symbol
        assert elemList[0].startswith(self._hdSyms[0]), "** ERROR: Sequence header doesn't start with '%s'" % self._hdSyms[0]
        #Validate the third line (quality header) strats with the expected symbol
        assert elemList[2].startswith(self._hdSyms[1]), "** ERROR: Quality header doesn't start with '%s'" % self._hdSyms[1]
        #Ensure that the length of the sequence (line2) matches the length of quality score (line 4)
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: Sequence and quality lengths do not match."
        #return the record as a tuple containing the four lines of FASTQ entry
        return tuple(elemList)
############### End of ParseFastQ Class ###############

def read_sample_barcodes(harrington_clinical_data):
    """Load sample barcodes from the clinical data file."""
    sample_barcodes = {}      #Initialize an empty dictionary to store sample names and barcodes
    with open(harrington_clinical_data, 'r') as file:  #open the clinical data file in read mode
        next(file)  # Skip the header line of the file
        for line in file:  #Iterate through each subsequent line in the file
            columns = line.strip().split()   # Remove whitespace/newlines and split the line into column
            sample_name = columns[0]   #Assume the first column contains the sample name
            barcode = columns[2]   # Assume the thired column contains the brcode
            sample_barcodes[sample_name] = barcode    #store the sample name and its barcode in the dictionary
    return sample_barcodes    # return the dictionary containing sample barcodes

def trim_barcode(seq, qual, barcode_length):
    """Trim the barcode and its associated quality scores from the read."""
    trimmed_seq = seq[barcode_length:]  #Remove the barcode  portion from the start of the sequence
    trimmed_qual = qual[barcode_length:]  # Remove the corresponding quality score
    return trimmed_seq, trimmed_qual   # Return the trimmed sequence and quality scores

def trim_low_quality(seq, qual):
    """Trim the read if two consecutive low-quality scores (D or F) are found at the end."""
    for i in range(len(qual) - 1):  # Loop through the quality scores until the second to last character
        if qual[i] in 'DF' and qual[i + 1] in 'DF':  #check if two consequtive low-quality scores are present
            return seq[:i], qual[:i] # Trim the sequence and quality scores up to this index
    return seq, qual   #If no low-quality tail is found , return the original sequence and quality scores

#  paths( Clinical dats, and output dictionary for demultiplexed FastQ files)
fastq_path = "./hawkins_pooled_sequences.fastq"  #  path to the pooled FastQ file containing all reads
clinical_file = "./harrington_clinical_data.txt"  # path to the clinical data file with sample barcode info
fastq_dir  = "./fastqs"  # output directory where demultiplexed FastQ will be stored

## create output directory
os.makedirs(fastq_dir, exist_ok=True)

# Load sample barcodes from clinical data file
sample_barcodes = read_sample_barcodes(clinical_file)
print("Loaded barcodes:", sample_barcodes)

# Create output file handles for each sample, one FASTQ file per sample
output_files = {sample: open(os.path.join(fastq_dir, f"{sample}.fastq"), "w") for sample in sample_barcodes.keys()}

# Parse the pooled FastQ file using the ParseFastQ class 
fastq_parser = ParseFastQ(fastq_path)
#Iterate over each FastQ 
for seq_header, seq, qual_header, qual in fastq_parser:
    #check the sample's barcode to determine which sample this read belongs to
    for sample, barcode in sample_barcodes.items():
        if seq.startswith(barcode):
            # Step 1: Trim the barcode from the sequence and quality strings
            trimmed_seq, trimmed_qual = trim_barcode(seq, qual, len(barcode))
            
            # Step 2: Trim low-quality ends from the read if applicable
            trimmed_seq, trimmed_qual = trim_low_quality(trimmed_seq, trimmed_qual)

            # Write the trimmed read to the appropriate sample's FastQ file
            output_files[sample].write(f"{seq_header}\n{trimmed_seq}\n{qual_header}\n{trimmed_qual}\n")
            break  # Move to the next read once a match is found

# Close all output files
for file in output_files.values():
    file.close()

##### End of Demultiplexing -section 1A and 1B


#### Start alignment on each FASTQ to the reference sequence ####
##### Start Alignment Pipeline #####

# Define paths for refrence genome and aklignment output directory
reference_genome = "./dgorgon_reference.fa"  # Reference genome
align_dir = "./bams"  # Directory for SAM, BAM, and sorted BAM

# Create output directory
os.makedirs(align_dir, exist_ok=True)

# Step 1: Index the reference genome (only if it hasn't been indexed)
if not os.path.exists(reference_genome + ".bwt"):  
    print(f"Indexing reference genome: {reference_genome}")
    subprocess.run(["bwa", "index", reference_genome], check=True)

# Step 2: Process each demultiplexed FASTQ file for alignment
for fastq_file in os.listdir(fastq_dir):
    if fastq_file.endswith(".fastq"):
        fastq_path = os.path.join(fastq_dir, fastq_file)
        sample_name = os.path.basename(fastq_file).replace(".fastq", "")
        sam_file = os.path.join(align_dir, f"{sample_name}.sam")
        bam_file = os.path.join(align_dir, f"{sample_name}.bam")
        sorted_bam_file = os.path.join(align_dir, f"{sample_name}.sorted.bam")

        # Align FASTQ to reference genome
        print(f"Aligning {fastq_file} to reference")
        with open(sam_file, "w") as sam_output:
            subprocess.run(["bwa", "mem", reference_genome, fastq_path], stdout=sam_output, check=True)

        # Convert SAM to BAM and Sort
        print(f"Converting {sam_file} to sorted BAM")
        subprocess.run(["samtools", "view", "-bS", sam_file, "-o", bam_file], check=True)
        subprocess.run(["samtools", "sort", "-o", sorted_bam_file, bam_file], check=True)

        # Index Sorted BAM
        print(f"Indexing {sorted_bam_file}")
        subprocess.run(["samtools", "index", sorted_bam_file], check=True)

        # Cleanup intermediate files
        print(f"Cleaning up {sam_file} and {bam_file}")
        os.remove(sam_file)
        os.remove(bam_file)

print("Alignment pipeline completed successfully.")

##### End of Alignment section

#### start  Variant Discovery

# Define paths
variant_dir = "./variant_results"  # Directory for storing variant discovery results

# Ensure variant results directory exists
os.makedirs(variant_dir, exist_ok=True)

def pileup_analysis(bam_file):
    """Perform pileup analysis on a sorted BAM file and detect variants."""
    sample_name = os.path.basename(bam_file).replace(".sorted.bam", "")
    variant_output_file = os.path.join(variant_dir, f"{sample_name}_variants.txt")

    # Open the BAM file 
    samfile = pysam.AlignmentFile(bam_file, "rb")
    ref_fasta = pysam.FastaFile(reference_genome)  ## Add refrence

    with open(variant_output_file, "w") as output:
        # Write header for output file
        output.write("Position\tCoverage\tReference\tBaseCounts\tBaseFrequencies\tVariantDetected\tMutation\n")

        # Since our reference has a single sequence, we process all reads
        for pileupcolumn in samfile.pileup():
            position = pileupcolumn.pos  # Reference position
            coverage = pileupcolumn.n  # Total reads covering this position

            ##retrieve The refrence base using the refrence name  from the BAM
            ref_name = samfile.get_reference_name(pileupcolumn.tid)
            ref_base = ref_fasta.fetch(ref_name, position, position+1)

            # Dictionary to count occurrences of each nucleotide
            ntdict = {}

            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # Uncomment the below line to debug the base calls
                    
                    print("\tBase in read %s = %s" % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))

                    base = pileupread.alignment.query_sequence[pileupread.query_position]

                    # Populate the dictionary with the base counts
                    ntdict[base] = ntdict.get(base, 0) + 1

            print(f"Coverage at base {position} = {coverage}")  # Debugging
            print(ntdict)  # Print base counts for debugging

            # Skip positions with no valid base calls
            if not ntdict:
                continue

            # Calculate base frequencies
            total_reads = sum(ntdict.values())
            base_frequencies = {nt: round((count / total_reads) * 100, 2) for nt, count in ntdict.items()}

            # Determine if there is a variant (more than one nucleotide present)
            variant_detected = len(base_frequencies) > 1


            ### Determine mutation call: if variant detected, select most common alternate allele
            if variant_detected:
                # Exclude the reference base when considering alternatives
                alt_bases = {nt: count for nt, count in ntdict.items() if nt != ref_base}
                if alt_bases:
                    alt_allele = max(alt_bases, key=alt_bases.get)
                    mutation = f"{ref_base}>{alt_allele}"
                else:
                    mutation = "No mutation"
            else:
                mutation = "No mutation"

            # Write results to file including the new columns
            output.write(f"{position}\t{coverage}\t{ref_base}\t{ntdict}\t{base_frequencies}\t{variant_detected}\t{mutation}\n")


           
    samfile.close()
    ref_fasta.close()
    print(f"Variant analysis completed for {bam_file}. Output saved to {variant_output_file}")

# Process each sorted BAM file for variant discovery
for sorted_bam_file in os.listdir(align_dir):
    if sorted_bam_file.endswith(".sorted.bam"):
        bam_path = os.path.join(align_dir, sorted_bam_file)

        # Run pileup analysis
        print(f"Performing variant discovery on {sorted_bam_file}")
        pileup_analysis(bam_path)

print("Variant discovery completed successfully.")




