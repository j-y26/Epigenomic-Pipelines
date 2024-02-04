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

**Table of contents**

- [Transcriptomics](#transcriptomics)
  - [Configuration](#configuration)
  - [Pre-alignment processing](#pre-alignment-processing)
    - [Quality control](#quality-control)
    - [Trimming](#trimming)
    - [Quality control after trimming](#quality-control-after-trimming)
  - [Alignment](#alignment)


## Configuration

Ensure that all the scripts in this folder are executable. If not, run the
following command:

```bash
chmod -R a+x ./transcriptome
```

Make sure to set up the proper parameters in the `config_rnaseq.sh` file. The
[configuration file](config_rnaseq.sh) is an example bash script that sets up
the parameters for the analysis. Most parameters are provided by default, but
the most important ones that must be set are:

- **`projPath`**: the path to the home project analysis folder
- **`fastqDir`**: the path to the folder containing the fastq files, should be a subfolder of `projPath`
- **`lane`**: the lane number of the sequencing run
- **`threads`**: the number of threads to use for the analysis

In addition, some program used are not linux, but Java-based, and require
setting the executable (`.jar` file) path in the `config_rnaseq.sh` file. 
For example, set `TRIMMOMATIC_JAR` to the path of the `trimmomatic-0.39.jar` file.

<br><br/>

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

The execution of all the scripts follow the same convention, where there should
be one and only one argument, which is the configuration file. Here, we run the
script as follows:

```bash
./fastqc_pre-trimming.sh config_rnaseq.sh
```

### Trimming

After checking the quality of the raw data, it is important to remove low quality
bases and adapter sequences. This can be done using tools such as Trimmomatic,
which is a flexible read trimming tool for Illumina NGS data.

Trimmomatic is a Java-based program, and the path to the `.jar` file must be set
in the `config_rnaseq.sh` file. In addition, the adapter sequences must be
provided to the program. While the program has a default set of adapter sequences,
located along with the `.jar` file, it is important to choose the correct
adapter sequences for the specific library preparation kit used. Consult the 
instruction for the kit to select the correct adapter sequences.

The trimming step make use of the `trimming_reads.sh` script, where we perform
trimming on the raw sequencing data.

```bash
./trimming_reads.sh config_rnaseq.sh
```

### Quality control after trimming

To assess whether the trimming has improved the quality of the reads, we perform
another round of quality control using FastQC. This is done using the
`fastqc_post-trimming.sh` script.

```bash
./fastqc_post-trimming.sh config_rnaseq.sh
```

We should expect to see an improvement in the quality of the reads after trimming.
This will facilitate a better alignment to the reference genome.

<br><br/>

## Alignment