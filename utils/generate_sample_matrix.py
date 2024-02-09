# Given a directory of files, find all files having the given suffix and
# generate a sample matrix csv file, with the following columns:
# - Sample
# - Label
# - Group
# - Replicate
# - Batch
# - Mark
# - PeakType
# - FileName

# Usage: python3 sample_matrix.py <directory> <suffix> <output_file>

import os
import sys
import pandas as pd

# Check if the number of arguments is correct
if len(sys.argv) != 4:
    print("Usage: python3 sample_matrix.py <directory> <suffix> <output_file>")
    sys.exit(1)

# Parse the command line arguments
directory = sys.argv[1]
suffix = sys.argv[2]
output_file = sys.argv[3]

# Traverse the directory to find files with the specified suffix
sample_matrix = []
for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith(suffix):
            # Extract sample name
            sample_name = file[:-len(suffix)]
            
            # Append sample data to the sample matrix
            sample_matrix.append({
                "Sample": sample_name,
                "Label": sample_name,
                "Group": "",
                "Replicate": "",
                "Batch": "",
                "Mark": "",
                "PeakType": "",
                "FileName": file
            })

# Convert the sample matrix to a pandas DataFrame
sample_matrix_df = pd.DataFrame(sample_matrix)
sample_matrix_df.to_csv(output_file, index=False)

print("Sample matrix generated successfully!")
print(f"Sample matrix saved to {output_file}")

# Exit the program
sys.exit(0)

# [END]