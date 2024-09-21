---
layout: default
title: Single-cell Multiome
nav_order: 5
---

# Single-cell Multiome

This page describes scripts that facilitate analyzing multiomic data at the
single cell resolution.

**Table of contents**

- [Single-cell Multiome](#single-cell-multiome)
  - [Gene Expression (GEX) Quality Control](#gene-expression-gex-quality-control)
    - [Proportion of exonic reads (exon\_prop)](#proportion-of-exonic-reads-exon_prop)


## Gene Expression (GEX) Quality Control

### Proportion of exonic reads (exon_prop)

The proportion of exonic reads for single nucleus sequencing is an important
matirx to assess the quality of the data. This is because that for GEM (Gel
bead in EMulsion) generation on single nucleus, only isolated nuclei are
captured in the GEMs. Therefore, RNA captured inside the nuclei are mostly
unspliced. High proportion of exonic reads indicates that the nuclei are
likely contaminated by cytoplasmic RNA due to broken nuclei.

Given the 10x Cell Ranger ARC alignment output (`gex_possorted_bam.bam`), the
tags inside the BAM files are used to indicate the properties of the reads.
The first script `gex_bam_tags_to_csv.py` extracts the tags from the BAM file
while only keeping the confidently mapped reads (MAPQ = 255) and have the read
collapsed to the unique molecular identifier (UMI).

```bash
python gex_bam_tags_to_csv.py <input_bam> <output_csv>
```

The second script `calc_read_type_prop.py` calculates the count and proportion
of reads that are exonic, intronic, and intergenic. The output is a CSV file
with the following columns:
- CB_cell_barcode: cell barcode
- total_reads: total number of confidently mapped reads
- exon_reads: number of confidently mapped exonic reads
- exon_prop: proportion of confidently mapped exonic reads
- intron_reads: number of confidently mapped intronic reads
- intron_prop: proportion of confidently mapped intronic reads
- intergenic_reads: number of confidently mapped intergenic reads
- intergenic_prop: proportion of confidently mapped intergenic reads

```bash
python calc_read_type_prop.py <input_csv> <output_csv> [chunk_size]
```

`chunk_size` is an optional argument that specifies the number of lines to read
from the input CSV file at a time, which is assigned to each worker process.
This script is parallelized using the `multiprocessing` module in Python.