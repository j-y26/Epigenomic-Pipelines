# Given a sample matrix csv file, extract the list of samples (labels) for each peak type.

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

# Print out the list of invalid samples
if not invalid_samples.empty:
    print('Invalid peak type annotation:')
    print(invalid_samples)
    print("Note that the peak type for samples must be either 'narrow' or 'broad'.")
    print("Samples annotated otherwise will not be used to call peaks.")
    print("An exception is negative control samples, which should not be annotated.")

# Keep only the samples with peak type 'narrow' or 'broad'
df = df[df['Peak_Type'].isin(['narrow', 'broad'])]

# Group the samples by peak type
grouped = df.groupby('Peak_Type')['Label'].apply(list)

# Write the list of samples for each peak type to a text file
for peak_type, samples in grouped.items():
    with open(os.path.join(sys.argv[2], f'samples_{peak_type}_peak.txt'), 'w') as f:
        for sample in samples:
            f.write(f'{sample}\t')

print("")
print('Samples extracted by peak type.')

# [END]