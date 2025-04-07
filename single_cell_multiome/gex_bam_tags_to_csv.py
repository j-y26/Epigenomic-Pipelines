# Parsing the cell ranger alignment tags from the bam files to csv files
# This script is used to parse the important cell ranger alignment tags from the
# gex_possorted_bam.bam file to a csv file.

# The output is a csv file such that:
#   - Each row corresponds to a read
#   - Each column corresponds to a tag
#   - The values are the tag values, without the tag name or type
#   - The UMIs are collapsed (i.e., no duplicates)
#   - Only confidently mapped reads are included (i.e., mapping quality = 255)
#   - Only cells with valid barcodes are included (CB and UB)

# The tags & types that are parsed are:
#   1. CB:Z: error-corrected and confirmed cell barcode
#   2. UB:Z: error-corrected and confirmed UMI
#   3. MAPQ (column 5): mapping quality
#   4. RE:A: single character indicating the region type of this alignment (E = exonic, N = intronic, I = intergenic)
#   5. GN:Z: gene name(s) for this alignment


# Usage: python gex_bam_tags_to_csv.py -i/--input <input_bam> -o/--output <output_csv> -t/--threads <num_threads>


# [Main]

import os
import pysam
import csv
import argparse
import multiprocessing as mp
from queue import Empty
from collections import defaultdict
import time

# Define BIN size for parallel processing
BIN_SIZE = 10000000  # 10 Mb bins

# Function to process a chunk of reads
def process_chunk(chunk):
    """
    Process a chunk of reads and return the parsed data and UMI dictionary.
    Each read is processed to extract the relevant tags, and the UMIs are collapsed.
    """
    umi_dict = defaultdict(set)
    results = []

    input_bam, reference, start_pos, end_pos = chunk
    
    # Open the BAM file and fetch reads from the chunk's region
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for read in bam.fetch(reference=reference, start=start_pos, end=end_pos):
            # Get the mapping quality, and if it is less than 255, skip the read
            mapq = read.mapping_quality
            if mapq < 255:
                continue

            # Get the cell barcode
            cb = read.get_tag("CB") if read.has_tag("CB") else ""
            if cb == "":
                continue

            # Get the UMI
            ub = read.get_tag("UB") if read.has_tag("UB") else ""
            if ub == "":
                continue

            # Add UMI to the dictionary
            umi_dict[cb].add(ub)

            # Get the region type
            re = read.get_tag("RE") if read.has_tag("RE") else ""

            # Write the tags to the results list
            gn = read.get_tag("GN") if read.has_tag("GN") else ""
            results.append([cb, ub, re, gn])

    return results, umi_dict


# Function to merge the UMI dictionaries from all threads
def merge_umi_dicts(umi_dicts):
    """
    Merge the UMI dictionaries from all threads into a single dictionary.
    The keys are the cell barcodes, and the values are the set of UMIs.
    """
    merged_dict = defaultdict(set)

    for d in umi_dicts:
        for cb, ub_set in d.items():
            merged_dict[cb].update(ub_set)
    return merged_dict


# Getting the chunks of the bam file
def get_bam_chunks(input_bam, num_threads, bin_size=BIN_SIZE):
    """
    Divide the BAM file into chunks for parallel processing.
    
    Returns:
        list: List of tuples (input_bam, reference, start_pos, end_pos)
    """
    chunks = []
    
    try:
        with pysam.AlignmentFile(input_bam, "rb") as bam:
            # Check if the BAM file is indexed
            try:
                if bam.has_index():
                    references = bam.references
                    references_lengths = bam.lengths
                    
                    # Create chunks based on reference sequences
                    for ref, length in zip(references, references_lengths):
                        for start in range(0, length, bin_size):
                            end = min(start + bin_size, length)
                            chunks.append((input_bam, ref, start, end))
                else:
                    raise ValueError("BAM file is not indexed")
            except (ValueError, AttributeError):
                # BAM is not indexed, process the whole file in chunks
                # Divide the file based on thread count
                print("BAM file is not indexed. Dividing file into regions based on thread count...")

                # Get approximate read count
                total_reads = 0
                sample_size = 1000
                sample_count = 0
                
                for read in bam.fetch(until_eof=True):
                    sample_count += 1
                    if sample_count > sample_size:
                        break
                
                # Estimate based on file size
                bam.reset()
                bam_filesize = os.path.getsize(input_bam)
                reads_per_thread = max(10000, bam_filesize // (10000 * num_threads))
                
                # Create chunks for the entire file
                for i in range(num_threads):
                    start_pos = i * reads_per_thread
                    end_pos = (i + 1) * reads_per_thread if i < num_threads - 1 else 0
                    chunks.append((input_bam, "", start_pos, end_pos))
                    
    except Exception as e:
        print(f"Error while analyzing BAM file: {e}")
        # Fallback to simple chunking
        chunks = [(input_bam, "", 0, 0)]
            
    return chunks


# The main function
def main(input_bam, output_csv, num_threads, bin_size=BIN_SIZE):
    """
    Main function to parse the cell ranger alignment tags from the bam file to a csv file.
    """
    start_time = time.time()

    # Check if the input bam file exists
    if not os.path.exists(input_bam):
        raise FileNotFoundError(f"Input BAM file '{input_bam}' does not exist.")
    
    print(f"Processing BAM file: {input_bam}")

    # Divide the BAM file into chunks for parallel processing
    chunks = get_bam_chunks(input_bam, BIN_SIZE)
    total_chunks = len(chunks)
    print(f"Divided BAM file into {total_chunks} chunks of size {BIN_SIZE} for parallel processing.")

    # Process the BAM file in parallel
    print(f"Starting {num_threads} threads for processing...")

    # Use a context manager for the pool to ensure proper cleanup
    with mp.Pool(processes=num_threads) as pool:
        results = pool.map(process_chunk, chunks)

    # Collect results
    all_results = []
    umi_dicts = []

    for result, umi_dict in results:
        all_results.extend(result)
        umi_dicts.append(umi_dict)
    
    print("Merging UMI dictionaries from all threads...")
    merged_umi_dict = merge_umi_dicts(umi_dicts)
    unique_cells = len(merged_umi_dict)
    total_umis = sum(len(umis) for umis in merged_umi_dict.values())
    
    print(f"Detected {unique_cells} cells with {total_umis} unique UMIs")
    
    print("Writing results to CSV file...")
    with open(output_csv, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["CB_cell_barcode", "UB_umi_barcode", "RE_region_type", "GN_gene_name"])
        writer.writerows(all_results)

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"Parsed {len(all_results)} reads successfully.")
    print(f"Finished writing results to {output_csv}.")
    print(f"Total processing time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="Parse the cell ranger alignment tags from the bam file to a csv file.")
    parser.add_argument("-i", "--input", required=True, help="Input bam file, must be sorted by coordinate")
    parser.add_argument("-o", "--output", required=True, help="Output csv file")
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count(), help="Number of threads to use")
    
    args = parser.parse_args()
    main(args.input, args.output, args.threads)

# [END]