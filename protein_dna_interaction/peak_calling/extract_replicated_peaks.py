# Given a sample matrix csv file, extract the list of samples (labels) that are
# biological replicates of each other. This allows the identification of
# replicates for each group of samples that are used to find consensus peaks.

# The sample csv file must have the following columns:
# - Label: The sample label
# - Peak_Type: The peak type of the sample

# Usage: python3 extract_samples_by_peak_type.py <sample_matrix> <output_dir>

import sys
import os
import pandas as pd

# Read in the sample matrix csv file
df = pd.read_csv(sys.argv[1])

# Make the Peak_Type column lowercase
df['Peak_Type'] = df['Peak_Type'].str.lower()

# Find all samples where Peak_Type is not narrow or broad
invalid_samples = df[~df['Peak_Type'].isin(['narrow', 'broad'])]

# Keep only the samples with peak type 'narrow' or 'broad'
df = df[df['Peak_Type'].isin(['narrow', 'broad'])]

# Group the samples by Group and Mark
grouped = df.groupby(['Group', 'Mark'])['Label'].apply(list)

# Write the list of samples for each group to a text file
for group, samples in grouped.items():
    with open(os.path.join(sys.argv[2], f'samples_{group[0]}_{group[1]}.txt'), 'w') as f:
        for sample in samples:
            f.write(f'{sample}\t')

print("")
print('Replicated samples extracted by group.')

# [END]