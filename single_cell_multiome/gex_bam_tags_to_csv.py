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
import tempfile
import hashlib
import time

BIN_SIZE = 10_000_000  # 10Mb bin for parallel processing

def process_chunk_to_temp(args):
    """
    Processes a chunk of a BAM file and writes results to a temporary CSV file.
    Returns path to temp file and a set of unique (CB, UB) hashes for global deduplication.
    """

    chunk, temp_dir = args
    input_bam, reference, start, end = chunk
    print(f"[Worker {mp.current_process().name}] Processing {reference}:{start}-{end}")

    os.makedirs(temp_dir, exist_ok=True)
    temp_fd, temp_path = tempfile.mkstemp(suffix=".tmp", prefix="chunk_", dir=temp_dir, text=True)
    os.close(temp_fd)

    umi_seen = set()

    with open(temp_path, "w", newline='') as f_out:
        writer = csv.writer(f_out)
        with pysam.AlignmentFile(input_bam, "rb") as bam:
            for read in bam.fetch(reference=reference, start=start, end=end):
                if read.mapping_quality != 255:
                    continue

                if not (read.has_tag("CB") and read.has_tag("UB")):
                    continue

                cb = read.get_tag("CB")
                ub = read.get_tag("UB")

                # Deduplicate within this chunk using a simple hash
                umi_key = f"{cb}:{ub}"
                umi_hash = hashlib.md5(umi_key.encode()).hexdigest()
                if umi_hash in umi_seen:
                    continue
                umi_seen.add(umi_hash)

                re = read.get_tag("RE") if read.has_tag("RE") else ""
                gn = read.get_tag("GN") if read.has_tag("GN") else ""

                writer.writerow([cb, ub, re, gn])

    return temp_path, umi_seen

def get_bam_chunks(input_bam, bin_size=BIN_SIZE):
    """
    Returns a list of (input_bam, reference, start, end) tuples for binning.
    """
    chunks = []
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        references = bam.references
        ref_lengths = bam.lengths
        for ref, length in zip(references, ref_lengths):
            for start in range(0, length, bin_size):
                end = min(start + bin_size, length)
                chunks.append((input_bam, ref, start, end))
    return chunks

def merge_temp_files(temp_paths, umi_hashes_global, final_output_path):
    """
    Merges all temp files and resolves UMI duplicates globally.
    """
    seen_global = set()
    with open(final_output_path, "w", newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(["CB_cell_barcode", "UB_umi_barcode", "RE_region_type", "GN_gene_name"])
        for temp_path in temp_paths:
            with open(temp_path, "r") as f:
                reader = csv.reader(f)
                for row in reader:
                    cb, ub = row[0], row[1]
                    key = f"{cb}:{ub}"
                    h = hashlib.md5(key.encode()).hexdigest()
                    if h not in seen_global:
                        seen_global.add(h)
                        writer.writerow(row)
            os.remove(temp_path)

def main(input_bam, output_csv, num_threads, temp_dir):
    start_time = time.time()
    print(f"Analysis started at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")

    if not os.path.exists(input_bam):
        raise FileNotFoundError(f"Input BAM file '{input_bam}' does not exist.")

    chunks = get_bam_chunks(input_bam)
    print(f"Divided BAM into {len(chunks)} chunks.")

    chunk_args = [(chunk, temp_dir) for chunk in chunks]

    with mp.Pool(processes=num_threads) as pool:
        results = pool.map(process_chunk_to_temp, chunk_args)

    temp_paths, umi_hashes_per_worker = zip(*results)
    total_umi_seen = set.union(*umi_hashes_per_worker)

    print(f"Merging {len(temp_paths)} chunk outputs...")
    print(f"Started at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))}")
    merge_temp_files(temp_paths, total_umi_seen, output_csv)

    print(f"Done at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))}")
    print(f"Unique UMIs: {len(total_umi_seen)}")
    print(f"Output written to: {output_csv}")

    print(f"Analysis finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))}")
    print(f"Total time: {time.time() - start_time:.2f}s")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Efficiently parse Cell Ranger BAM tags into a CSV.")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file (must be coordinate-sorted and indexed).")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file.")
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count(), help="Number of threads to use.")
    parser.add_argument("-d", "--temp_dir", default=tempfile.gettempdir(), help="Temporary directory for intermediate files.")
    args = parser.parse_args()
    main(args.input, args.output, args.threads, args.temp_dir)

# [END]