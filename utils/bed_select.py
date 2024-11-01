# Given an input BED file with multiple sequence coordinates, extract and output 
# a subset of sequences based on regular expression patterns matching a list of 
# sequence identifiers in user-specified columns to the BED file.

# Usage: python bed_select.py -i/--input <input_bed_file>
#                             -o/--output <output_bed_file>
#                             -c/--columns <column1> [<column2> ...]
#                             -p/--pattern [<pattern>]
#                             -pf/--pattern-file [<pattern_file>]
#                             -e/--exclude [<exclude_pattern>]
#                             -ef/--exclude-file [<exclude_pattern_file>]

import argparse
import re
import pandas as pd

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

def select_sequences(input, output, in_patterns, ex_patterns, columns):
    """
    Selects sequences based on the inclusion and exclusion patterns and writes the output to a new BED file.
    """
    print("Begin selecting sequences...")

    # Read the input BED file
    bed_df = pd.read_csv(input, sep="\t", header=None)

    # Convert 1-based columns to 0-based for indexing
    columns = [col - 1 for col in columns]

    # Vectorized inclusion filtering: check if any pattern matches in specified columns
    if in_patterns:
        inclusion_mask = pd.concat([bed_df[col].astype(str).str.contains(pat, regex=True) for col in columns for pat in in_patterns], axis=1).any(axis=1)
    else:
        inclusion_mask = pd.Series(True, index=bed_df.index)

    # Vectorized exclusion filtering: check if any pattern matches in specified columns
    if ex_patterns:
        exclusion_mask = pd.concat([bed_df[col].astype(str).str.contains(pat, regex=True) for col in columns for pat in ex_patterns], axis=1).any(axis=1)
    else:
        exclusion_mask = pd.Series(False, index=bed_df.index)
    
    # Apply the inclusion and exclusion masks to select the sequences
    selected_df = bed_df[inclusion_mask & ~exclusion_mask]
    print(f"Selected {selected_df.shape[0]} sequences.")

    # Print the selected sequences at the specified columns
    print(selected_df[columns])

    # Check if there is any sequence selected
    if selected_df.empty:
        print("No matching records found. Output BED file not created.")
    else:
        # Write the selected sequences to the output BED file
        selected_df.to_csv(output, sep="\t", header=False, index=False)
        print(f"Selected sequences written to {output}")

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="Given an input BED file with multiple sequence coordinates, extract and output a subset of sequences based on regular expression patterns matching a list of sequence identifiers in user-specified columns to the BED file.")

    parser.add_argument("-i", "--input", required=True, help="Input BED file, no header allowed")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument("-c", "--columns", nargs="+", type=int, default=[4], \
                        help="Column(s) containing sequence identifiers, 1-based index, must be a positive integer. Default: 4")
    parser.add_argument("-p", "--pattern", dest="pattern", help="Pattern to include")
    parser.add_argument("-pf", "--pattern-file", dest="pattern_file", help="File containing patterns to include")
    parser.add_argument("-e", "--exclude", dest="exclude_pattern", help="Pattern to exclude")
    parser.add_argument("-ef", "--exclude-file", dest="exclude_pattern_file", help="File containing patterns to exclude")

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

    # Identify the column(s) in the BED file to search for patterns
    print("Columns to search for patterns:")
    print(args.columns)

    # Select sequences based on the inclusion and exclusion patterns
    select_sequences(args.input, args.output, in_patterns, ex_patterns, args.columns)

if __name__ == "__main__":
    main()