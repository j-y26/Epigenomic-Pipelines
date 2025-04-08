# Calculating the proportion of reads of each type in each cell barcode

# This script is used subsequently to the gex_bam_tags_to_csv.py script by
# operating on the output csv file.

# This script runs on multi-core CPUs to speed up the processing of the data
# by dividing the data into chunks and processing each chunk in parallel.

# Read proportion are calculated as:
#     exon_prop = number of confidently mapped exonic reads / number of confidently mapped reads
#     intron_prop = number of confidently mapped intronic reads / number of confidently mapped reads
#     intergenic_prop = number of confidently mapped intergenic reads / number of confidently mapped reads

# Assumptions of the input csv file:
#   - The input csv file has the following columns:
#       - CB_cell_barcode: cell barcode
#       - UB_umi_barcode: UMI barcode
#       - RE_region_type: region type of the alignment
#       - GN_gene_name: gene name(s) for the alignment
#   - The input csv file only consists of confidently mapped reads (MAPQ >= 255)
#   - Only reads with valid cell barcodes (CB) and UMI barcodes (UB) are included
#   - The UMIs are collapsed (i.e., no duplicates)

# The output is a csv file such that:
#   - Each row corresponds to a cell
#   - Each column corresponds to metadata, with columns
#       - CB_cell_barcode: cell barcode
#       - total_reads: total number of confidently mapped reads
#       - exon_reads: number of confidently mapped exonic reads
#       - exon_prop: proportion of confidently mapped exonic reads
#       - intron_reads: number of confidently mapped intronic reads
#       - intron_prop: proportion of confidently mapped intronic reads
#       - intergenic_reads: number of confidently mapped intergenic reads
#       - intergenic_prop: proportion of confidently mapped intergenic reads

# Usage: python calc_read_type_prop.py -i <input_csv> -o <output_csv> -c [<chunk_size>]

# [Main]

import sys
import os
import csv
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import itertools
import argparse

# Define default dictionary structure
def default_dict_int():
    return {"E": 0, "N": 0, "I": 0}

# Handle read type counts for each chunk
def process_chunk(chunk):
    # Exonic (E), Intronic (N), Intergenic (I)
    print(f"Processing new chunk of {len(chunk)} reads... on process {os.getpid()}")
    read_type_dict = defaultdict(default_dict_int)
    for row in chunk:
        try:
            cell_barcode, _, region_type, _ = row
            if region_type in ["E", "N", "I"]:
                read_type_dict[cell_barcode][region_type] += 1
        except Exception as e:
            print(f"Error processing row: {row}. Error: {str(e)}")
    return dict(read_type_dict)

# chuck iterator
def chunk_iterator(reader, chunk_size):
    while True:
        chunk = list(itertools.islice(reader, chunk_size))
        if not chunk:
            break
        yield chunk

# Main function
def main(input_csv, output_csv, chunk_size):
    # Check if the input csv file exists
    if not os.path.exists(input_csv):
        print("The input csv file does not exist.")
        sys.exit(1)

    print("Begin reading input csv file...")

    combined_dict = defaultdict(default_dict_int)
    
    # Open the input csv file for reading
    with open(input_csv, "r") as f:
        reader = csv.reader(f)
        header = next(reader) # Skip the header
        
        # Assign chucks to be processed by available CPUs while reading the file
        available_cpus = cpu_count() - 1 
        if available_cpus < 1:
            available_cpus = 1
        
        print(f"Reading and processing the input csv file with {available_cpus} CPUs on the fly...")
        
        with Pool(processes=available_cpus) as pool:
            for result in pool.imap_unordered(process_chunk, chunk_iterator(reader, chunk_size)):
                for cell_barcode, read_type_counts in result.items():
                    for region_type, count in read_type_counts.items():
                        combined_dict[cell_barcode][region_type] += count
    

    # Print the number of reads
    total_cell_barcodes = len(combined_dict)
    print(f"Completed processing {total_cell_barcodes} unique cell barcodes.")
    
    # Calculate the total read counts and proportions
    # Run calculations and writing on the fly
    with open(output_csv, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["CB_cell_barcode", "total_reads", "exon_reads", "exon_prop", "intron_reads", "intron_prop", "intergenic_reads", "intergenic_prop"])
        
        for cell_barcode, read_type_counts in combined_dict.items():
            total_read_count = sum(read_type_counts.values())
            exon_read_count = read_type_counts["E"]
            intron_read_count = read_type_counts["N"]
            intergenic_read_count = read_type_counts["I"]

            exon_prop = exon_read_count / total_read_count if total_read_count > 0 else 0
            intron_prop = intron_read_count / total_read_count if total_read_count > 0 else 0
            intergenic_prop = intergenic_read_count / total_read_count if total_read_count > 0 else 0

            writer.writerow([cell_barcode, total_read_count, 
                             exon_read_count, exon_prop, 
                             intron_read_count, intron_prop, 
                             intergenic_read_count, intergenic_prop])

    print("Done.")
    print(f"Output written to {output_csv}")

# Check if the script is being run directly
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Calculate the proportion of reads of each type in each cell barcode.")
    parser.add_argument("-i", "--input_csv", required=True, help="Input csv file")
    parser.add_argument("-o", "--output_csv", required=True, help="Output csv file")
    parser.add_argument("-c", "--chunk_size", type=int, default=1000000, help="Chunk size for processing (default: 1000000)")
    args = parser.parse_args()

    # Run the main function with the provided arguments
    main(args.input_csv, args.output_csv, args.chunk_size)

# [END]
