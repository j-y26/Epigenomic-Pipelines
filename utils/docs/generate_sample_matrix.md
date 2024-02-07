---
layout: default
title: Sample Matrix
parent: Utilities
nav_order: 3
---

# Generate sample matrix

The `generate_sample_matrix.py` script is used to generate a sample matrix from 
the file names in a given directory. In Illumina sequencing, it is common
to use file names to identify samples, but very often that the files names
are not easy to understand, and may not reflect the properties of the samples.

Here, given a directory containing files, the script will generate a sample matrix
with the following columns:
- `Sample`: the original sample name
- `Label`: the user-defined label
- `Group`: the user-defined group
- `Replicate`: the user-defined replicate
- `Batch`: the user-defined batch
- `Mark`: the user-defined mark
- `PeakType`: the user-defined peak type
- `FileName`: the original file name

We will use all files in the directory matching the given file suffix.

## Usage

```bash
python3 sample_matrix.py <directory> <suffix> <output_file>
```

For example, given a directory `./data` containing the following files:
  
```
sample1_rep1_trimmed.fastq.gz
sample2_rep2_trimmed.fastq.gz
sample3_rep1_trimmed.fastq.gz
```

We can generate a sample matrix with the following command:

```bash
python3 sample_matrix.py ./data _trimmed.fastq.gz sample_matrix.csv
```

The generated `sample_matrix.csv` will look like this:

| Sample | Label | Group | Replicate | Batch | Mark | PeakType | FileName |
|--------|-------|-------|-----------|-------|------|----------|----------|
| sample1_rep1 | | | | | | | sample1_rep1_trimmed.fastq.gz |
| sample2_rep2 | | | | | | | sample2_rep2_trimmed.fastq.gz |
| sample3_rep1 | | | | | | | sample3_rep1_trimmed.fastq.gz |

This sample matrix is applicable to most type of analysis. However, some columns
may not be applicable to some analysis, and the user can choose to remove them
from the sample matrix.