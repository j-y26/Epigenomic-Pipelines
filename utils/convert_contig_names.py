# Convert chromosome names in a input file based on a mapping file

# Usage: python3 convert_contig_names.py <input_file> <mapping_file> <output_file> [<chromosome_column>]

import sys
import pandas as pd

input_file = sys.argv[1]
mapping_file = sys.argv[2]
output_file = sys.argv[3]
if len(sys.argv) > 4:
    chromosome_column = int(sys.argv[4]) - 1
else:
    chromosome_column = 0

# Read mapping file
mapping = pd.read_csv(mapping_file, sep='\t', header=None, dtype={0: str, 1: str})
mapping.columns = ['old', 'new']

# Read input file
input = pd.read_csv(input_file, sep='\t', header=None, dtype={chromosome_column: str})

# Convert chromosome names
input[chromosome_column] = input[chromosome_column].map(mapping.set_index('old')['new'])

# Write output
input.to_csv(output_file, sep='\t', header=False, index=False)

print('Done!')

# [END]