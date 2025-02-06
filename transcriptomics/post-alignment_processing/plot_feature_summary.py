# Summarize feature summary results and plot the summary metrics.

# Usage: python plot_feature_summary.py -i/--input <feature_count_summary_file>
#                                       -s/--suffix [<bam_file_suffix>]
#                                       -w/--width [<plot_width>]
#                                       -ht/--height [<plot_height>]


import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend before importing pyplot
import matplotlib.pyplot as plt

# Function to plot alignment metrics
def plot_metrics(input_file, bam_suffix, plot_width, plot_height):
    """Plot feature summary metrics."""
    # Read the featureCount summary file
    # Check if the input file exists
    if not os.path.exists(input_file):
        sys.exit('Error: File not found: ' + input_file)
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    # Sample names in header are paths to the BAM files
    # Remove path and suffix to get sample names
    df.columns = df.columns.str.split('/').str[-1]  # Remove path
    df.columns = df.columns.str.replace(bam_suffix, '')  # Remove suffix

    # Use millions as unit
    df = df / 1e6

    # Transpose the dataframe
    df = df.T

    # Calculate total analyzed reads
    df['Total_reads'] = df.sum(axis=1)

    # Reset names
    df = df.reset_index().rename(columns={'index': 'Sample'})

    # Plot metrics
    fig, axes = plt.subplots(3, 2, figsize=(plot_width, plot_height))
    fig.suptitle('FeatureCounts summary metrics')

    # Define color function
    def get_colors(condition, true_color='#d62728', false_color='#1f77b4'):
        """Return a list of colors where condition is met; otherwise, use default."""
        return [true_color if cond else false_color for cond in condition]
    
    # Define bbox for text with transparent background and border
    def bbox_props(color, alpha=0.7):
        return dict(facecolor=color, alpha=alpha, edgecolor='grey', boxstyle="round,pad=0.3")
    
    # Plot feature summary metrics
    df.plot(x='Sample', y='Assigned', kind='bar', ax=axes[0, 0], title='Assigned reads (Millions)',
            color=get_colors(df['Assigned'] < 0.5 * df['Total_reads']))  # Red if < 50% of total reads
    axes[0, 0].set_xlabel("")
    axes[0, 0].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[0, 0].axhline(y=df['Assigned'].mean(), color='k', linestyle='--', label='Mean')
    axes[0, 0].text(0.02, 0.07, "Highlighted if < 50% of total reads", ha='left', va='top', transform=axes[0, 0].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Total_reads', kind='bar', ax=axes[0, 1], title='Total alignments (Millions)',
            color=get_colors((df['Total_reads'] < 0.6 * df['Total_reads'].mean()) | (df['Total_reads'] > 1.4 * df['Total_reads'].mean())))  # Red if < 60% or > 140% of the mean
    axes[0, 1].set_xlabel("")
    axes[0, 1].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[0, 1].axhline(y=df['Total_reads'].mean(), color='k', linestyle='--', label='Mean')
    axes[0, 1].text(0.02, 0.07, "Highlighted if < 60% or > 140% of the mean", ha='left', va='top', transform=axes[0, 1].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Unassigned_MultiMapping', kind='bar', ax=axes[1, 0], title='Unassigned: Multi-mapping (Millions)',
            color=get_colors(df['Unassigned_MultiMapping'] > 0.3 * df['Total_reads']))  # Red if > 30% of total reads
    axes[1, 0].set_xlabel("")
    axes[1, 0].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[1, 0].axhline(y=df['Unassigned_MultiMapping'].mean(), color='k', linestyle='--', label='Mean')
    axes[1, 0].text(0.02, 0.07, "Highlighted if > 30% of total reads", ha='left', va='top', transform=axes[1, 0].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Unassigned_Ambiguity', kind='bar', ax=axes[1, 1], title='Unassigned: Ambiguity (Millions)',
            color=get_colors(df['Unassigned_Ambiguity'] > 0.15 * df['Total_reads']))  # Red if > 15% of total reads
    axes[1, 1].set_xlabel("")
    axes[1, 1].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[1, 1].axhline(y=df['Unassigned_Ambiguity'].mean(), color='k', linestyle='--', label='Mean')
    axes[1, 1].text(0.02, 0.07, "Highlighted if > 15% of total reads", ha='left', va='top', transform=axes[1, 1].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Unassigned_NoFeatures', kind='bar', ax=axes[2, 0], title='Unassigned: No features (Millions)',
            color=get_colors(df['Unassigned_NoFeatures'] > 0.05 * df['Total_reads']))  # Red if > 5% of total reads
    axes[2, 0].set_xlabel("")
    axes[2, 0].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[2, 0].axhline(y=df['Unassigned_NoFeatures'].mean(), color='k', linestyle='--', label='Mean')
    axes[2, 0].text(0.02, 0.07, "Highlighted if > 5% of total reads", ha='left', va='top', transform=axes[2, 0].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Unassigned_Secondary', kind='bar', ax=axes[2, 1], title='Unassigned: Secondary (Millions)',
            color=get_colors(df['Unassigned_Secondary'] > 0.05 * df['Total_reads'])) # Red if > 5% of total reads
    axes[2, 1].set_xlabel("")
    axes[2, 1].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[2, 1].axhline(y=df['Unassigned_Secondary'].mean(), color='k', linestyle='--', label='Mean')
    axes[2, 1].text(0.02, 0.07, "Highlighted if > 5% of total reads", ha='left', va='top', transform=axes[2, 1].transAxes, fontsize=9, bbox=bbox_props('white'))

    plt.tight_layout()
    plt.savefig(input_file + '.png')
    plt.close()
    

# Main
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Summarize STAR alignment results')
    parser.add_argument('-i', '--input', help='Path to the feature count summary file', required=True)
    parser.add_argument('-s', '--suffix', help='BAM file suffix, default: _name_sorted.bam', default='_name_sorted.bam')
    parser.add_argument('-w', '--width', help='Plot width, default: 10', default=10)
    parser.add_argument('-ht', '--height', help='Plot height, default: 12', default=12)
    args = parser.parse_args()

    # Read the featureCount summary file
    plot_metrics(args.input, args.suffix, args.width, args.height)

    print('Feature summary metrics plot saved to', args.input + '.png')

if __name__ == '__main__':
    main()