# Given a csv file with the following columns: Sample, Label
# sample: name of the sample
# label: label of the sample
# Replace the "sample" part of the file name with the label

# Usage: python file_name_conversion.py <csv_file> <file_dir> <file_extension>

# Note: the file names should follow the format: sample<file_extension>
# For example: sample_bowtie2.bam has the extension _bowtie2.bam

import os
import sys
import pandas as pd

csv_file = sys.argv[1]
file_dir = sys.argv[2]
file_extension = sys.argv[3]

# Read in the csv file
sample_sheet = pd.read_csv(csv_file)

# Extract the sample names and labels   
sample_names = sample_sheet["Sample"]
sample_labels = sample_sheet["Label"]

# Rename the files, by iterating through all the files in the directory
for file in os.listdir(file_dir):
    # Check if the file has the correct extension
    if file.endswith(file_extension):
        # Obtain the sample name, by removing the extension
        sample_name = file.replace(file_extension, "")
        # Check if the sample name is in the sample sheet
        if sample_name in sample_names.values:
            # Obtain the index of the sample name
            sample_index = sample_names[sample_names == sample_name].index[0]
            # Obtain the label of the sample
            sample_label = sample_labels[sample_index]
            # Rename the file
            os.rename(os.path.join(file_dir, file), os.path.join(file_dir, sample_label + file_extension))
            print("Renamed " + file + " to " + sample_label + file_extension)
        else:
            print("Sample " + sample_name + " not found in the sample sheet")
    else:
        print("File " + file + " does not have the correct extension")

# [END]