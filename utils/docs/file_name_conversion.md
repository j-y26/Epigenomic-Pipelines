---
layout: default
title: Convert file names
parent: Utilities
nav_order: 2
---

# Batch file name conversion

In many cases, sample names used during sequencing might not be easy for
identifying samples during the analysis, particularly if user want to label
samples by their group properties (e.g., treatment, time point, etc.), or 
biological replicates. Here, we provide a utility script to replace the sample
names in the file names with user-defined label.

Given a comma-separated (csv) file containing sample names (used in the files)
and a user-refined label, the script will replace the sample names in the file
names while remaining the file extension for all files in a given directory.

It is important to note that the script must contain two columns:
- `Sample`: the original sample name, should match the sample names in the file names
- `Label`: the user-defined label to replace the sample name

Using this csv file, users can regard it has a sample tracking sheet, and
can contain additional metadata information for each sample. For example, in
ChIP-seq experiments, the csv file can contain the following columns:
- `Sample`: the original sample name
- `Label`: the user-defined label
- `Treatment`: the treatment condition
- `Batch`: the batch number
- `Replicate`: the biological replicate number
- `Mark`: the antibody used

This is performed using the `file_name_conversion.py` script, which can be found
[here](https://github.com/j-y26/Epigenomic-Pipelines/blob/main/utils/file_name_conversion.py).

To use the script, users can run the following command:

```bash
python3 file_name_conversion.py <path_to_csv_file> <path_to_directory> <file_extension>
```

For example: `s1_H3K27_bowtie2.bam` could be a bad file naming for analysis.
Given a csv file containing the following information:

| Sample | Label | Treatment | Batch | Replicate | Mark |
|--------|-------|-----------|-------|-----------|------|
| s1_H3K27 | ctrl1_H3K27me3_rep1 | Control | 1 | 1 | H3K27me3 |

Running the script with the following command:

```bash
python3 file_name_conversion.py sample_tracking.csv /path/to/directory/ _bowtie2.bam
```

will rename the file to `ctrl1_H3K27me3_rep1_bowtie2.bam`.