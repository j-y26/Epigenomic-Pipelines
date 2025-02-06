# Summarize STAR alignment results and generate a CSV file with a plot of the alignment metrics.

# Usage: python summarize_alignment.py -d/--directory [<star_output_directory>]
#                                      -s/--suffix [<alignment_log_file_suffix>]
#                                      -o/--output [<output_csv_file>]
#                                      -w/--width [<plot_width>]
#                                      -ht/--height [<plot_height>]

import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend before importing pyplot
import matplotlib.pyplot as plt


# Function to parse STAR alignment log file
def parse_star_log(log_file, suffix):
    """Parse STAR alignment log file and return a dictionary."""
    metrics = {
        'Sample': os.path.basename(log_file).replace(suffix, ''),
        'Total_reads': None,
        'Average_input_read_length': None,
        'Uniquely_mapped_reads_number': None,
        'Uniquely_mapped_reads_percent': None,
        'Number_of_reads_mapped_to_multiple_loci': None,
        'Percent_of_reads_mapped_to_multiple_loci': None,
        'Number_of_reads_mapped_to_too_many_loci': None,
        'Percent_of_reads_mapped_to_too_many_loci': None,
        'Number_of_reads_unmapped_too_many_mismatches': None,
        'Percent_of_reads_unmapped_too_many_mismatches': None,
        'Number_of_reads_unmapped_too_short': None,
        'Percent_of_reads_unmapped_too_short': None,
        'Number_of_reads_unmapped_other': None,
        'Percent_of_reads_unmapped_other': None,
        'Average_mapped_length': None,
        'Mismatch_rate_per_base': None,
        'Deletion_rate_per_base': None,
        'Insertion_rate_per_base': None,
    }

    with open(log_file, 'r') as f:
        for line in f:
            parts = line.strip().split('|')
            if len(parts) < 2:
                continue
            key, value = parts[0].strip(), parts[1].strip()

            if key == 'Number of input reads':
                metrics['Total_reads'] = int(value)
            elif key == 'Average input read length':
                metrics['Average_input_read_length'] = int(value)
            elif key == 'Uniquely mapped reads number':
                metrics['Uniquely_mapped_reads_number'] = int(value)
            elif key == 'Uniquely mapped reads %':
                metrics['Uniquely_mapped_reads_percent'] = float(value.replace('%', ''))
            elif key == 'Number of reads mapped to multiple loci':
                metrics['Number_of_reads_mapped_to_multiple_loci'] = int(value)
            elif key == '% of reads mapped to multiple loci':
                metrics['Percent_of_reads_mapped_to_multiple_loci'] = float(value.replace('%', ''))
            elif key == 'Number of reads mapped to too many loci':
                metrics['Number_of_reads_mapped_to_too_many_loci'] = int(value)
            elif key == '% of reads mapped to too many loci':
                metrics['Percent_of_reads_mapped_to_too_many_loci'] = float(value.replace('%', ''))
            elif key == 'Number of reads unmapped: too many mismatches':
                metrics['Number_of_reads_unmapped_too_many_mismatches'] = int(value)
            elif key == '% of reads unmapped: too many mismatches':
                metrics['Percent_of_reads_unmapped_too_many_mismatches'] = float(value.replace('%', ''))
            elif key == 'Number of reads unmapped: too short':
                metrics['Number_of_reads_unmapped_too_short'] = int(value)
            elif key == '% of reads unmapped: too short':
                metrics['Percent_of_reads_unmapped_too_short'] = float(value.replace('%', ''))
            elif key == 'Number of reads unmapped: other':
                metrics['Number_of_reads_unmapped_other'] = int(value)
            elif key == '% of reads unmapped: other':
                metrics['Percent_of_reads_unmapped_other'] = float(value.replace('%', ''))
            elif key == 'Average mapped length':
                metrics['Average_mapped_length'] = float(value)
            elif key == 'Mismatch rate per base, %':
                metrics['Mismatch_rate_per_base'] = float(value.replace('%', ''))
            elif key == 'Deletion rate per base':
                metrics['Deletion_rate_per_base'] = float(value.replace('%', ''))
            elif key == 'Insertion rate per base':
                metrics['Insertion_rate_per_base'] = float(value.replace('%', ''))

    return metrics

# Function to plot alignment metrics
def plot_metrics(df, plt_name, plot_width, plot_height):
    """Plot alignment metrics from STAR log files."""

    # Calculate total mapped reads
    df['Total_reads_M'] = df['Total_reads'] / 1e6
    df['Total_mapped_percent'] = df['Uniquely_mapped_reads_percent'] + df['Percent_of_reads_mapped_to_multiple_loci']
    df['Total_mapped_reads'] = df['Total_reads_M'] * df['Total_mapped_percent'] / 100

    
    fig, axes = plt.subplots(3, 2, figsize=(plot_width, plot_height))
    fig.suptitle('STAR Alignment metrics')

    # Define color function
    def get_colors(condition, true_color='#d62728', false_color='#1f77b4'):
        """Return a list of colors where condition is met; otherwise, use default."""
        return [true_color if cond else false_color for cond in condition]  # None allows default colors
    
    # Define bbox for text with transparent background and border
    def bbox_props(color, alpha=0.7):
        return dict(facecolor=color, alpha=alpha, edgecolor='grey', boxstyle="round,pad=0.3")
    
    # Plot alignment metrics
    df.plot(x='Sample', y='Total_reads_M', kind='bar', ax=axes[0, 0], title='Total input reads (Millions)', 
            color=get_colors(df['Total_reads'] < 0.6 * df['Total_reads'].mean()))  # Red if < 60% of mean
    axes[0, 0].set_xlabel("")
    axes[0, 0].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[0, 0].axhline(y=df['Total_reads_M'].mean(), color='k', linestyle='--', label='Mean')
    axes[0, 0].text(0.02, 0.07, "Highlighted if < 60% of mean sequencing depth", ha='left', va='top', transform=axes[0, 0].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Total_mapped_reads', kind='bar', ax=axes[0, 1], title='Total mapped reads (Millions)',
            color=get_colors(df['Total_mapped_percent'] < 80))  # Red if < 80%
    axes[0, 1].set_xlabel("")
    axes[0, 1].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[0, 1].axhline(y=df['Total_mapped_reads'].mean(), color='k', linestyle='--', label='Mean')
    axes[0, 1].text(0.02, 0.07, "Highlighted if < 80% of total reads are mapped", ha='left', va='top', transform=axes[0, 1].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Uniquely_mapped_reads_percent', kind='bar', ax=axes[1, 0], title='Uniquely mapped reads (%)',
            color=get_colors(df['Uniquely_mapped_reads_percent'] < 70))  # Red if < 70%
    axes[1, 0].set_xlabel("")
    axes[1, 0].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[1, 0].axhline(y=df['Uniquely_mapped_reads_percent'].mean(), color='k', linestyle='--', label='Mean')
    axes[1, 0].text(0.02, 0.07, "Highlighted if < 70% of reads are uniquely mapped", ha='left', va='top', transform=axes[1, 0].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Percent_of_reads_mapped_to_multiple_loci', kind='bar', ax=axes[1, 1], title='Mapped to multiple loci (%)',
            color=get_colors(df['Percent_of_reads_mapped_to_multiple_loci'] > 30))  # Red if > 30%
    axes[1, 1].set_xlabel("")
    axes[1, 1].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[1, 1].axhline(y=df['Percent_of_reads_mapped_to_multiple_loci'].mean(), color='k', linestyle='--', label='Mean')
    axes[1, 1].text(0.02, 0.07, "Highlighted if > 30% of reads are mapped to multiple loci", ha='left', va='top', transform=axes[1, 1].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Percent_of_reads_mapped_to_too_many_loci', kind='bar', ax=axes[2, 0], title='Mapped to too many loci (%)',
            color=get_colors(df['Percent_of_reads_mapped_to_too_many_loci'] > 1))  # Red if > 1
    axes[2, 0].set_xlabel("")
    axes[2, 0].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[2, 0].axhline(y=df['Percent_of_reads_mapped_to_too_many_loci'].mean(), color='k', linestyle='--', label='Mean')
    axes[2, 0].text(0.02, 0.07, "Highlighted if > 1% of reads are mapped to too many loci", ha='left', va='top', transform=axes[2, 0].transAxes, fontsize=9, bbox=bbox_props('white'))

    df.plot(x='Sample', y='Percent_of_reads_unmapped_too_short', kind='bar', ax=axes[2, 1], title='Unmapped: too short (%)',
            color=get_colors(df['Percent_of_reads_unmapped_too_short'] > 5))  # Red if > 5%
    axes[2, 1].set_xlabel("")
    axes[2, 1].get_legend().set_visible(False)
    # Add a horizontal line at the mean
    axes[2, 1].axhline(y=df['Percent_of_reads_unmapped_too_short'].mean(), color='k', linestyle='--', label='Mean')
    axes[2, 1].text(0.02, 0.07, "Highlighted if > 5% of reads are unmapped: too short", ha='left', va='top', transform=axes[2, 1].transAxes, fontsize=9, bbox=bbox_props('white'))  

    plt.tight_layout()
    plt.savefig(plt_name)
    plt.close()

# Function to summarize alignment results
def summarize_alignment(star_output_directory, alignment_log_file_suffix, output_file, plot_width, plot_height):
    """Summarizes STAR alignment results from log files into a CSV."""
    alignment_logs = [f for f in os.listdir(star_output_directory) if f.endswith(alignment_log_file_suffix)]
    
    if not alignment_logs:
        sys.exit(f'No alignment log files found with suffix {alignment_log_file_suffix} in {star_output_directory}')
    
    summary_data = []
    for log in alignment_logs:
        log_path = os.path.join(star_output_directory, log)
        summary_data.append(parse_star_log(log_path, alignment_log_file_suffix))
    
    df = pd.DataFrame(summary_data)
    df.to_csv(output_file, index=False)

    # Plot alignment metrics
    if output_file.endswith('.csv'):
        plt_name = output_file.replace('.csv', '.png')
    else:
        plt_name = output_file + '.png'
    
    plot_metrics(df, plt_name, plot_width, plot_height)

    print(f'Summary saved to {output_file}')
    print(f'Alignment metrics plot saved to {plt_name}')


# Main
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Summarize STAR alignment results')
    parser.add_argument('-d', '--directory', help='STAR output directory, default: ./star_output', default='./star_output')
    parser.add_argument('-s', '--suffix', help='BAM file suffix, default: Log.final.out', default='Log.final.out')
    parser.add_argument('-o', '--output', help='Output csv file, default: alignment_summary.csv', default='alignment_summary.csv')
    parser.add_argument('-w', '--width', help='Plot width, default: 10', default=10)
    parser.add_argument('-ht', '--height', help='Plot height, default: 12', default=12)
    args = parser.parse_args()

    # Summarize alignment results
    summarize_alignment(args.directory, args.suffix, args.output, args.width, args.height)

if __name__ == '__main__':
    main()