---
layout: default
title: Transcriptomics
nav_order: 2
---

# Transcriptomics

This page describes the steps to analyze RNA-seq data. We assume that the data
are generated following a standard RNA-seq protocol, with libraries sequenced
using the Illumina platform.

Here, we describe the analysis from the input fastq files.

## Configuration

Ensure that all the scripts in this folder are executable. If not, run the
following command:

```bash
chmod a+x -R transcriptome
```

Make sure to set up the proper parameters in the `config_rnaseq.sh` file. The
[configuration file](config_rnaseq.sh) is an example bash script that sets up
the parameters for the analysis. Most parameters are provided by default, but
the most important ones that must be set are:

- **`projPath`**: the path to the home project analysis folder
- **`fastqDir`**: the path to the folder containing the fastq files, should be a subfolder of `projPath`
- **`lane`**: the lane number of the sequencing run
- **`threads`**: the number of threads to use for the analysis

## Pre-alignment processing

### Quality control

Before starting the analysis, it is important to check the quality of the raw
sequencing data. This can be done using tools such as FastQC, which provides
information about the quality of the reads, including per base sequence quality,
per base sequence content, sequence length distribution, and overrepresented
sequences.

Although it is necessary to examine the sequence duplication levels, it is
important to note that for downstream differential expression analysis,
duplications are usually **preserved**, and the removal of duplicates is not
recommended.

The first QC step make use of the `fastqc_pre-trimming.sh` script, where we
perform QC on the raw sequencing data.