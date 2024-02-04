# Extracting gene information from gtf file
# Input: gtf file
# Output: gene information file (gene_id, gene_name, gene_type, chromosome, gene_start, gene_end, gene_strand), tab separated
# Usage: python extract_gene_from_gtf.py <gtf_file> <output_file>

import sys

gtf_file = sys.argv[1]
output_file = sys.argv[2]

with open(gtf_file, 'r') as f:
    with open(output_file, 'w') as o:
        # write header
        o.write('gene_id' + '\t' + 'gene_name' + '\t' + 'gene_type' + '\t' + 'chromosome' + '\t' + 'gene_start' + '\t' + 'gene_end' + '\t' + 'gene_strand' + '\n')

        # Iterate through each line in the gtf file
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    gene_id = line[8].split(';')[0].split(' ')[1].replace('"', '')
                    gene_name = line[8].split(';')[2].split(' ')[2].replace('"', '')
                    gene_type = line[8].split(';')[4].split(' ')[2].replace('"', '')
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
                    gene_start = line[3]
                    gene_end = line[4]
                    gene_strand = line[6]
                    o.write(gene_id + '\t' + gene_name + '\t' + gene_type + '\t' + chromosome + '\t' + gene_start + '\t' + gene_end + '\t' + gene_strand + '\n')