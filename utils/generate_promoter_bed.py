# Given a transcript annotation file (generated by utils/extract_transcript_from_gtf.py)
# find the promoter regions of each transcript with user defined regions upstream and downstream of the TSS
# Output a bed file with the promoter regions, with gene id and gene symbol as the name field

# Usage: python generate_promoter_bed.py -t <transcript_file> -o <output_file> -u <upstream> -d <downstream>

import argparse
import pandas as pd
import numpy as np

def generate_promoter_bed(transcript_file, output_file, upstream, downstream):
    # Read in the transcript annotation file
    transcripts = pd.read_csv(transcript_file, sep="\t", header=0)

    # Keep only protein coding genes
    transcripts = transcripts[transcripts["transcript_type"] == "protein_coding"]

    # Extract the transcript_id, gene_name, chromosome, strand, start and end positions
    transcripts = transcripts[["transcript_id", "gene_name", "chromosome", "transcript_strand", "transcript_start", "transcript_end"]]

    # Generate the promoter regions
    conditions = [transcripts["transcript_strand"] == "+", transcripts["transcript_strand"] == "-"]
    choices_start = [transcripts["transcript_start"] - upstream, transcripts["transcript_end"] - downstream]
    choices_end = [transcripts["transcript_start"] + downstream, transcripts["transcript_end"] + upstream]

    transcripts["promoter_start"] = np.select(conditions, choices_start)
    transcripts["promoter_end"] = np.select(conditions, choices_end)

    # If any of the start or end positions are negative, set them to 0
    transcripts["promoter_start"] = transcripts["promoter_start"].clip(lower=0)
    transcripts["promoter_end"] = transcripts["promoter_end"].clip(lower=0)

    # Generate the bed file
    bed = transcripts[["chromosome", "promoter_start", "promoter_end", "gene_name", "transcript_id", "transcript_strand"]]

    # Write the bed file
    bed.to_csv(output_file, sep="\t", header=False, index=False)

def main():
    parser = argparse.ArgumentParser(
        description='Generate promoter regions from a transcript annotation file\n Please use utils/extract_transcript_from_gtf.py to generate the transcript annotation file'
    )
    
    parser.add_argument(
        '-t', '--transcript_file', 
        help='Transcript annotation file, generated by running utils/extract_transcript_from_gtf.py on a gtf file', 
        required=True)
    
    parser.add_argument(
        '-o', '--output_file', 
        help='Output bed file', 
        required=True)
    
    parser.add_argument(
        '-u', '--upstream', 
        type=int, 
        default=3000, 
        help='Number of bases upstream of the TSS')
    
    parser.add_argument(
        '-d', '--downstream', 
        type=int, default=1000, 
        help='Number of bases downstream of the TSS')
    
    parser.print_help()
    args = parser.parse_args()

    # Generate the promoter regions
    generate_promoter_bed(args.transcript_file, args.output_file, args.upstream, args.downstream)

    print("Promoter regions BED file generated: " + args.output_file)

if __name__ == "__main__":
    main()