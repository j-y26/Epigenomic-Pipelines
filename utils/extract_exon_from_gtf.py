# Extracting exon information from gtf file
# Input: gtf file
# Output: exon information file (exon_id, exon_number, exon_start, exon_end, exon_strand, chromosome, gene_id, gene_name), tab separated
# Usage: python extract_exon_from_gtf.py <gtf_file> <output_file>

import sys
import re

gtf_file = sys.argv[1]
output_file = sys.argv[2]

with open(gtf_file, 'r') as f:
    with open(output_file, 'w') as o:
        # write header
        o.write('exon_id' + '\t' + 'exon_number' + '\t' + 'exon_start' + '\t' + 'exon_end' + '\t' + 'strand' + '\t' + 'chromosome' + '\t' + 'gene_id' + '\t' + 'gene_name' + '\n')

        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if line[2] == 'exon':
                    # Extract exon_id using regular expression
                    exon_id_match = re.search(r'exon_id "([^"]+)"', line[8])
                    if exon_id_match:
                        exon_id = exon_id_match.group(1)

                    # Extract exon_number using regular expression
                    exon_number_match = re.search(r'exon_number "([^"]+)"', line[8])
                    if exon_number_match:
                        exon_number = exon_number_match.group(1)

                    chromosome = line[0]
                    # if chromosome is a number, add 'chr'
                    # if chromosome is a == 'X' or 'Y', add 'chr'
                    # if chromosome is a 'MT', add 'chr'
                    if chromosome == 'MT':
                        chromosome = 'chrM'
                    elif chromosome == 'X' or chromosome == 'Y':
                        chromosome = 'chr' + chromosome
                    elif chromosome.isdigit():
                        chromosome = 'chr' + chromosome
                    exon_start = line[3]
                    exon_end = line[4]
                    exon_strand = line[6]

                    # Gene information
                    gene_id = line[8].split(';')[0].split(' ')[1].replace('"', '')
                    gene_name_match = re.search(r'gene_name "([^"]+)"', line[8])
                    if gene_name_match:
                        gene_name = gene_name_match.group(1)

                    # Write to output file
                    o.write(exon_id + '\t' + exon_number + '\t' + exon_start + '\t' + exon_end + '\t' + exon_strand + '\t' + chromosome + '\t' + gene_id + '\t' + gene_name + '\n')
