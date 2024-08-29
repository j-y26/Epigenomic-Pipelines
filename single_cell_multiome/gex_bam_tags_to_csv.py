# Parsing the cell ranger alignment tags from the bam files to csv files
# This script is used to parse the important cell ranger alignment tags from the
# gex_possorted_bam.bam file to a csv file.

# The output is a csv file such that:
#   - Each row corresponds to a read
#   - Each column corresponds to a tag
#   - The values are the tag values, without the tag name or type
#   - The UMIs are collapsed (i.e., no duplicates)

# The tags & types that are parsed are:
#   1. CB:Z: error-corrected and confirmed cell barcode
#   2. UB:Z: error-corrected and confirmed UMI
#   3. MAPQ (column 5): mapping quality
#   4. RE:A: single character indicating the region type of this alignment (E = exonic, N = intronic, I = intergenic)
#   5. GN:Z: gene name(s) for this alignment


# Usage: python gex_bam_tags_to_csv.py <input_bam> <output_csv>


# [Main]

import sys
import os
import pysam
import csv

# Check if the input arguments are correct
if len(sys.argv) != 3:
    print("Usage: python gex_bam_tags_to_csv.py <input_bam> <output_csv>")
    sys.exit(1)

# Input arguments
input_bam = sys.argv[1]
output_csv = sys.argv[2]

# Check if the input bam file exists
if not os.path.exists(input_bam):
    print("The input bam file does not exist.")
    sys.exit(1)

# Open the input bam file for reading
bam = pysam.AlignmentFile(input_bam, "rb")

# Temporary dictionary to store the UMIs for each cell barcode
# The keys are the cell barcodes, and the values are the set of UMIs, which are used to collapse the UMIs
umi_dict = {}

# Open the output csv file for writing
with open(output_csv, "w") as f:
    writer = csv.writer(f)

    # Write the header
    writer.writerow(["CB_cell_barcode", "UB_umi_barcode", "MAPQ", "RE_region_type", "GN_gene_name"])

    # Iterate through the reads
    for read in bam.fetch():
        # Get the cell barcode
        cb = read.get_tag("CB") if read.has_tag("CB") else ""

        # Get the UMI
        ub = read.get_tag("UB") if read.has_tag("UB") else ""

        # Check if the cell barcode is already in the dictionary
        if cb in umi_dict:
            # Check if the UMI is already in the set
            if ub in umi_dict[cb]:
                continue
            else:
                umi_dict[cb].add(ub)
        else:
            umi_dict[cb] = set([ub])

        # Get the mapping quality
        mapq = read.mapping_quality

        # Get the region type
        re = read.get_tag("RE") if read.has_tag("RE") else ""

        # Write the tags to the csv file
        gn = read.get_tag("GN") if read.has_tag("GN") else ""

        writer.writerow([cb, ub, mapq, re, gn])

# Close the bam file
bam.close()

print("Finished parsing the tags from the bam file.")
print("Output csv file:", output_csv)

# [END]