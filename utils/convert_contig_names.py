# Convert chromosome names in a bed file based on a mapping file

# Usage: python3 convert_contig_names.py <bed_file> <mapping_file> <output_file>

import sys
import pandas as pd

bed_file = sys.argv[1]
mapping_file = sys.argv[2]
output_file = sys.argv[3]

# Read mapping file
mapping = pd.read_csv(mapping_file, sep='\t', header=None)
mapping.columns = ['old', 'new']

# Read bed file
bed = pd.read_csv(bed_file, sep='\t', header=None)

# Convert chromosome names
bed[0] = bed[0].map(mapping.set_index('old')['new'])

# Write output
bed.to_csv(output_file, sep='\t', header=False, index=False)

print('Done!')

# [END]