# Given a csv file for sample metadata, for each unique mark,
# generate a text file containing the sample labels for that mark.

# Note, the sample matrix must contain a "Label" and "Mark" column.

# Usage: python3 samples_for_mark.py <sample_matrix.csv> <output_dir>

import pandas as pd
import sys
import os

# Check that the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python3 samples_for_mark.py <sample_matrix.csv> <output_dir>")
    sys.exit(1)

# Read in the sample matrix
sample_matrix = pd.read_csv(sys.argv[1])

# Get a list of unique marks
marks = sample_matrix['Mark'].unique()

# Create a directory to store the output files
output_dir = sys.argv[2]
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# For each unique mark, generate a text file containing the sample labels for that mark
for mark in marks:
    mark_samples = sample_matrix[sample_matrix['Mark'] == mark]['Label']
    mark_samples_file = os.path.join(output_dir, f"samples_{mark}.txt")
    with open(mark_samples_file, 'w') as f:
        f.write('\t'.join(mark_samples))
    print(f"Generated {mark_samples_file}")

# [END]