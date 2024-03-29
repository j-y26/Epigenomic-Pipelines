# Given a csv file generated with mt_reads.sh that contains the number of
# mitochondrial and total reads for each sample, this script will generate a
# a bar plot to visualize the percentage of mitochondrial reads.
# A plot is generated for each mark.

# Usage: python3 mt_read_plot.py <csv file> <output_dir>

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# Read in the csv file
df = pd.read_csv(sys.argv[1])

# Calculate the percentage of mitochondrial reads
df['mito_pct'] = df['mt_reads'] / df['total_reads'] * 100

# Iterate through each mark and plot the percentage of mitochondrial reads
for mark in df['mark'].unique():
    mark_df = df[df['mark'] == mark]
    plt.bar(mark_df['sample'], mark_df['mito_pct'])
    plt.xlabel('Sample')
    plt.ylabel('Mitochondrial Reads (%)')
    plt.title(mark)

    # Rotate the x-axis labels
    plt.xticks(rotation=45)

    # Label the y-axis with percentages
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:.0f}%".format(x)))

    # Automatically adjust the layout
    plt.tight_layout()

    # Save the plot to pdf
    plt.savefig(os.path.join(sys.argv[2], f'mt_read_plot_{mark}.pdf'))

    # Close the plot
    plt.close()