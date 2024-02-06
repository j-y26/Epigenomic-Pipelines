---
layout: default
title: Chromatin Accessibility
nav_order: 4
---

# Chromatin Accessibility (ATAC-seq)

Chromatin accessibility is a measure of the ability of nuclear proteins to
interact with DNA. It is a key factor in the regulation of gene expression and
is often used to identify regulatory elements in the genome. Here, we describe
the processing and analysis of raw data files for ATAC-seq and DNase-seq
experiments.

**Table of contents**

- [Chromatin Accessibility (ATAC-seq)](#chromatin-accessibility-atac-seq)
  - [Configuration](#configuration)
  - [Pre-alignment processing](#pre-alignment-processing)
    - [Quality control](#quality-control)
    - [Trimming](#trimming)
    - [Quality control after trimming](#quality-control-after-trimming)

## Configuration

Ensure that all the scripts in this folder are executable. If not, run the
following command:

```bash
chmod -R a+x ./chromatin_accessibility
```

Make sure to set up the proper parameters in the `config_atacseq.sh` file. The
[configuration file](config_atacseq.sh) is an example bash script that sets up
the parameters for the analysis. Most parameters are provided by default, but
the most important ones that must be set are:

- **`projPath`**: the path to the home project analysis folder
- **`fastqDir`**: the path to the folder containing the fastq files, should be a subfolder of `projPath`
- **`lane`**: the lane number of the sequencing run
- **`threads`**: the number of threads to use for the analysis

In addition, some program used are not linux, but Java-based, and require
setting the executable (`.jar` file) path in the `config_chipseq.sh` file.
For example, set `PICARD_JAR` to the path of the `picard.jar` file.

<br><br/>

## Pre-alignment processing

The pre-alignment processing of ChIP-seq, CUT&Tag, and CUT&RUN data includes the
following steps:

- [FastQC processing of pre-trimmed reads](pre-alignment_processing/fastqc_pre-trimming.sh)
- [Trimming of reads with `fastp`](pre-alignment_processing/trimming_reads.sh)
- [FastQC processing of post-trimmed reads](pre-alignment_processing/fastqc_post-trimming.sh)

Here, we perform a step by step analysis to ensure the quality of the raw data
and to remove any adapter contamination or low-quality reads.

### Quality control

In the first step, we perform quality control of the raw reads using `fastqc`.
It is important to note that in genomic analysis, sequence duplication levels are
critical to assess. A lower sequence duplication level is preferred, and a high
duplication rate means a reduced complexity of the library. It is better to
optimize the final PCR amplification step to avoid bias and associated cost
due to the over-amplification of the library. In later steps, we will remove the
duplicates using `picard` tools so that only unique reads are used for the
alignment and the enrichment of "peaks".

The execution of all the scripts follow the same convention, where there should
be one and only one argument, which is the configuration file. Here, we run the
script as follows:

```bash
./fastqc_pre-trimming.sh config_atacseq.sh
```

### Trimming

After checking the quality of the raw data, it is important to remove low quality
bases and adapter sequences. This can be done using tools such as `fastp`, which
detects and cuts adapters, poly-N, and low-quality sequences.

```bash
./trimming_reads.sh config_atacseq.sh
```

### Quality control after trimming

After trimming the reads, we can again assess the quality of the trimmed reads
and examine whether the trimming procedure has improved the quality of the reads.
We should expect low quality bases are removed, along with adapter sequences.

```bash
./fastqc_post-trimming.sh config_atacseq.sh
```
