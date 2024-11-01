# Given an input FASTA file with multiple sequences, extract and output a subset
# of sequences based on regular expression patterns matching a list of sequence
# identifiers to the fasta header lines.

# Usage: python multi_fasta_select.py -i/--input <input_fasta_file>
#                                     -o/--output <output_fasta_file>
#                                     -p/--pattern [<pattern>]
#                                     -pf/--pattern-file [<pattern_file>]
#                                     -e/--exclude [<exclude_pattern>]
#                                     -ef/--exclude-file [<exclude_pattern_file>]

import argparse
import re
from Bio import SeqIO

def parse_patterns(pattern, pattern_file):
    """
    Parses the pattern(s) provided by the user and returns a list of compiled regex patterns.
    """
    patterns = []
    if pattern:
        print(pattern)
        patterns.append(re.compile(pattern))
    if pattern_file:
        with open(pattern_file, "r") as f:
            for line in f:
                print(line.strip())
                patterns.append(re.compile(line.strip()))
    return patterns

def select_sequences(input_fasta, output_fasta, in_patterns, ex_patterns):
    """
    Selects sequences from the input FASTA file based on inclusion and exclusion patterns
    and writes them to the output FASTA file.
    """
    print("Begin selecting sequences...")

    # Temporary variables to store matched patterns
    found_matching_records = False
    matched_records = []

    # Iterate through the input FASTA file and select sequences based on the inclusion and exclusion patterns
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.id

        # Check if the header matches any of the inclusion patterns
        include = any(pat.search(header) for pat in in_patterns) if in_patterns else True

        # Check if the header matches any of the exclusion patterns
        exclude = any(pat.search(header) for pat in ex_patterns) if ex_patterns else False

        # Write the sequence to the output FASTA file if it matches the inclusion pattern and does not match the exclusion pattern
        if include and not exclude:
            matched_records.append(record)
            found_matching_records = True
            print(f"Selected: {header}")
    print("FASTA sequence selection complete.")

    # Write the selected sequences to the output FASTA file
    if found_matching_records:
        with open(output_fasta, "w") as out_fasta:
            SeqIO.write(matched_records, out_fasta, "fasta")
        print(f"Selected sequences written to {output_fasta}")
    else:
        print("No sequences matched the given patterns. No output file generated.")

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="Given an input FASTA file with multiple sequences, extract and output a subset of sequences based on regular expression patterns matching a list of sequence identifiers to the FASTA header lines.")
    parser.add_argument("-i", "--input", dest="input_fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", dest="output_fasta", required=True, help="Output FASTA file")
    parser.add_argument("-p", "--pattern", dest="pattern", help="Pattern to match in the sequence identifier. Must provide either -p/--pattern or -pf/--pattern-file", default=None)
    parser.add_argument("-pf", "--pattern-file", dest="pattern_file", help="File containing patterns to match in the sequence identifier, where each regex pattern should begin on a new line. Overwrites -p/--pattern if provided", default=None)
    parser.add_argument("-e", "--exclude", dest="exclude_pattern", help="[Optional] Pattern to exclude from the sequence identifier", default=None)
    parser.add_argument("-ef", "--exclude-file", dest="exclude_pattern_file", help="[Optional]File containing patterns to exclude from the sequence identifier, where each regex pattern should begin on a new line. Overwrites -e/--exclude if provided", default=None)
    args = parser.parse_args()

    # Check if at least one of "pattern", "pattern_file" is provided
    if not args.pattern and not args.pattern_file:
        parser.error("At least one of --pattern or --pattern_file must be provided")
    
    # Check if both "pattern" and "pattern_file" are provided,
    # and if so, "pattern" will be ignored
    if args.pattern and args.pattern_file:
        print("Both --pattern and --pattern_file provided. Ignoring --pattern.")
        args.pattern = None

    # Check if both "exclude" and "exclude_file" are provided,
    # and if so, "exclude" will be ignored
    if args.exclude_pattern and args.exclude_pattern_file:
        print("Both --exclude and --exclude_file provided. Ignoring --exclude.")
        args.exclude_pattern = None

    # Parse the inclusion and exclusion patterns
    print("Parsing patterns...")
    print("Patterns to include:")
    in_patterns = parse_patterns(args.pattern, args.pattern_file)
    print("Patterns to exclude:")
    ex_patterns = parse_patterns(args.exclude_pattern, args.exclude_pattern_file)
    select_sequences(args.input_fasta, args.output_fasta, in_patterns, ex_patterns)

if __name__ == "__main__":
    main()