# Extracting transcript information from gtf file
# Input: gtf file
# Output: transcript information file (transcript_id, transcript_name, transcript_type, chromosome, transcript_start, transcript_end, transcript_strand, gene_id, gene_name), tab separated
# Usage: python extract_transcript_from_gtf.py <gtf_file> <output_file>

import sys
import re

gtf_file = sys.argv[1]
output_file = sys.argv[2]

with open(gtf_file, 'r') as f:
    with open(output_file, 'w') as o:
        # write header
        o.write('transcript_id' + '\t' + 'transcript_name' + '\t' + 'transcript_type' + '\t' + 'chromosome' + '\t' + 'transcript_start' + '\t' + 'transcript_end' + '\t' + 'transcript_strand' + '\t' + 'gene_id' + '\t' + 'gene_name' + '\n')

        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if line[2] == 'transcript':
                    # Extract transcript_id using regular expression
                    transcript_id_match = re.search(r'transcript_id "([^"]+)"', line[8])
                    if transcript_id_match:
                        transcript_id = transcript_id_match.group(1)

                    # Extract transcript_name using regular expression
                    transcript_name_match = re.search(r'transcript_name "([^"]+)"', line[8])
                    if transcript_name_match:
                        transcript_name = transcript_name_match.group(1)
                    
                    # Extract transcript_type using regular expression
                    transcript_type_match = re.search(r'transcript_biotype "([^"]+)"', line[8])
                    if transcript_type_match:
                        transcript_type = transcript_type_match.group(1)

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
                    transcript_start = line[3]
                    transcript_end = line[4]
                    transcript_strand = line[6]

                    # Gene information
                    gene_id = line[8].split(';')[0].split(' ')[1].replace('"', '')
                    gene_name_match = re.search(r'gene_name "([^"]+)"', line[8])
                    if gene_name_match:
                        gene_name = gene_name_match.group(1)

                    # Write to output file
                    o.write(transcript_id + '\t' + transcript_name + '\t' + transcript_type + '\t' + chromosome + '\t' + transcript_start + '\t' + transcript_end + '\t' + transcript_strand + '\t' + gene_id + '\t' + gene_name + '\n')

