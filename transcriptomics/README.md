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
  - [Running the entire pipeline](#running-the-entire-pipeline)
  - [Pre-alignment processing](#pre-alignment-processing)
    - [ORA decompression](#ora-decompression)
    - [Quality control](#quality-control)
    - [Trimming](#trimming)
    - [Quality control after trimming](#quality-control-after-trimming)
    - [rRNA Contamination Analysis and Cleaning \[Optional\]](#rrna-contamination-analysis-and-cleaning-optional)
    - [Quality control after rRNA cleaning](#quality-control-after-rrna-cleaning)
  - [Alignment](#alignment)
    - [Building STAR index](#building-star-index)
    - [Alignment](#alignment-1)
  - [Post-alignment processing](#post-alignment-processing)
    - [Sorting BAM file by name](#sorting-bam-file-by-name)
    - [Generating count matrix](#generating-count-matrix)
    - [Visualizing the RNA-seq Tracks](#visualizing-the-rna-seq-tracks)


## Configuration

Ensure that all the scripts in this folder are executable. If not, run the
following command:

```bash
chmod -R a+x ./transcriptomics
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

## Running the entire pipeline

We understand that a standardized RNA-seq processing pipeline is important to
remain consistent for all samples in different batches/groups. In the meantime,
in high-performance computing (HPC) environments, constant submission of jobs
can be cumbersome. Therefore, we provide a script that runs the entire pipeline
in one go.

This involves a sequential execution of the scripts in the following order:

1. Quality control on the raw sequencing data
2. Trimming the reads
3. Quality control after trimming
4. Aligning the reads to the reference genome
5. Sorting the BAM file by name
6. Generating the count matrix
7. Converting the BAM file to bigWig format

To run the entire pipeline, we make use of the `&&` operator, which requires
subsequent commands to run only if the previous command was successful. This
ensures that the pipeline is run sequentially.

```bash
~/pre-alignment_processing/fastqc_pre-trimming.sh config_rnaseq.sh && \
~/pre-alignment_processing/trimming_reads_fastp.sh config_rnaseq.sh && \
~/pre-alignment_processing/fastqc_post-trimming.sh config_rnaseq.sh && \
~/alignment/star_alignment.sh config_rnaseq.sh && \
~/post-alignment_processing/sort_bam_by_name.sh config_rnaseq.sh && \
~/post-alignment_processing/featureCounts_gene.sh config_rnaseq.sh && \
~/post-alignment_processing/bam_to_coverage.sh config_rnaseq.sh
```

Note:

- The `&&` operator ensures that the subsequent command is run only if the
  previous command was successful. Hence, always check the output of each
  command to ensure that it ran successfully.
- The `config_rnaseq.sh` file must be set up correctly before running the 
  pipeline. The paths to the reference genome index file and GTF annotation
  file must be set in the configuration file.
- This script assumes that a valid STAR index has been built for the reference
  genome. If not, run the `build_star_index.sh` script before running the
  pipeline.

We then describe the details of each step in the pipeline.

<br><br/>

## Pre-alignment processing

### ORA decompression

Some sequencing centers provide the raw sequencing data in the ORA format, which
provides a higher compression ratio compared to the standard gzip format. The
ORA format can be decompressed using the `ora_to_gz.sh` script.

```bash
./ora_to_gz.sh config_rnaseq.sh
```

The script utilizes the `orad` program provided by Illumina. For details,
refer to the [Illumina website](https://support.illumina.com/sequencing/sequencing_software/DRAGENORA.html).
The `DRAGEN ORA` suite requires a reference to decompress the ORA files. The
reference can be downloaded from the [Illumina website reference file releases](https://support.illumina.com/downloads/ora-decompression-reference-files.html).

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

The trimming step make use of the `trimming_reads_fastp.sh` script, where we perform
trimming on the raw sequencing data.

```bash
./trimming_reads_fastp.sh config_rnaseq.sh
```

Trimmomatic is a Java-based program, and the path to the `.jar` file must be set
in the `config_rnaseq.sh` file. In addition, the adapter sequences must be
provided to the program. While the program has a default set of adapter sequences,
located along with the `.jar` file, it is important to choose the correct
adapter sequences for the specific library preparation kit used. Consult the 
instruction for the kit to select the correct adapter sequences.

Here, we use fastp because it is faster than Trimmomatic and provides similar
performance. In addition, the program provides polyG and polyX trimming, which
is useful for removing artifacts in the reads. PolyG can arise in the NextSeq
and NovaSeq platforms due to a lack of signal at the ends of the reads (read2).
PolyT can arise due to the polyA tail in the mRNA or issue with the oligo-dT
beads during library preparation.

### Quality control after trimming

To assess whether the trimming has improved the quality of the reads, we perform
another round of quality control using FastQC. This is done using the
`fastqc_post-trimming.sh` script.

```bash
./fastqc_post-trimming.sh config_rnaseq.sh
```

We should expect to see an improvement in the quality of the reads after trimming.
This will facilitate a better alignment to the reference genome.

### rRNA Contamination Analysis and Cleaning [Optional]

Ribosomal RNA (rRNA) contamination is a common issue in RNA-seq data, especially
when the RNA is not polyA-selected, i.e., in rRNA-depleted samples. rRNA
depletion efficiency can vary, especially in low-quality RNA samples, and
contamination can affect the downstream analysis, particularly the estimation
of gene expression levels.

To check for rRNA contamination, we use `SortMeRNA`, which is a program that
can be used to filter rRNA reads from the RNA-seq data. The program requires
a reference database of rRNA sequences, which can be downloaded from the
`SortMeRNA` website. Furthermore, we also use species-specific annotation of
rRNA sequences, which can be downloaded from databases such as Ensembl or UCSC.
To extract the `FASTA` sequences of the rRNA genes, we use first select the
rRNA genes from the GTF annotation file and then extract the sequences from
the reference genome FASTA file. This is easily achieved using the `grep` and
`bedtools getfasta` commands.

To check for rRNA contamination, we use the `sortmerna_rRNA.sh` script.

```bash
./sortmerna_align.sh config_rnaseq.sh
```

The script will output aligned and unaligned reads. The program also outputs the
proportion of reads that align to the rRNA database. Unaligned reads are in
the `fastq` format and can be used for downstream analysis, which are cleaned
reads that do not contain rRNA contamination.

An alignment file is also generated in the BAM format, which can be visualized
using tools such as `IGV` to check the alignment of the reads to the rRNA genes.

### Quality control after rRNA cleaning

After cleaning the rRNA reads, it is important to check the quality of the
cleaned reads. This can be done using tools such as FastQC, which provides
information about the quality of the reads, including per base sequence quality,
per base sequence content, sequence length distribution, and overrepresented
sequences. A shift of GC content towards the center should be expected after
rRNA cleaning, since rRNA reads are usually GC-rich.

The quality control step after rRNA cleaning make use of the `fastqc_post-sortmerna.sh`
script.

```bash
./fastqc_post-sortmerna.sh config_rnaseq.sh
```

<br><br/>

## Alignment

Aligning the reads to the reference genome is the next step in the analysis. This
can be done using tools such as `STAR`, which is a fast and accurate aligner for
RNA-seq data. However, `STAR` is memory-intensive, and the amount of memory
required depends on the size of the reference genome. For example, the human
reference genome requires approximately 30GB of memory.

### Building STAR index

Before running the alignment, it is important to prepare the reference genome
index. This can be done using the `STAR --runMode genomeGenerate` command.
You can download the reference genome FASTA file and GTF annotation file from
a database such as Ensembl or UCSC. Alternatively, STAR provides a set of
pre-built reference genome indexes for many species, which can be downloaded
from the STAR website.

To build the index of the reference genome, we use the `build_star_index.sh`
script.

```bash
./build_star_index.sh config_rnaseq.sh
```

Ensure that you have the path to the reference genome FASTA file and GTF
annotation file set in the `config_rnaseq.sh` file. The `STAR` program will
use these files to build the index and output the index files to the
directory you specified.

### Alignment

After building the index, we can now align the reads to the reference genome.
Here we use the `star_alignment.sh` script.

```bash
./star_alignment.sh config_rnaseq.sh
```

The script will align the reads to the reference genome and output the aligned
reads in the BAM format. The BAM file is a binary format, and can be converted
to a human-readable SAM format using the `samtools view` command. Furthermore,
the default behavior of the pipeline ensures that the BAM file is sorted by
coordinate, which are usually required for indexing the BAM file and for
some downstream analysis.

In case where the alignment fails to the having very short fragments 
(Umapped: too short), this is primarily due to the fact that RNA degradation is
severe and the RNA fragments are too short to be aligned. In such cases, it is
recommended to use the `--outFilterMatchNminOverLread` parameter to increase the
minimum match length. This parameter specifies the minimum fraction of the read
length that must match the reference genome. The default value is 0.66, which
means that at least 66% of the read length must match the reference genome. If
the RNA fragments are very short, this parameter can be decreased.

The `star_alignment_short_fragment.sh` script can be used to align the reads
with short fragments. However, parameters need to be optimized.

```bash
./star_alignment_short_fragment.sh config_rnaseq.sh
```

<br><br/>

## Post-alignment processing

After aligning the reads to the reference genome, it is important to check the
quality of the alignment. In particular, `STAR` provides a set of quality
metrics that can be used to assess the quality of the alignment, particularly,
the alignment rate, the number of uniquely mapped reads, and the number of
multi-mapped reads.

Furthermore, we need to translate the alignment files into a count matrix,
which can be used for downstream differential expression analysis.

### Sorting BAM file by name

To use `featureCounts` to generate the count matrix, it is important to sort the
BAM file by name. This can be done using the `samtools sort` command.

The `sort_bam_by_name.sh` script sorts the BAM file by name.

```bash
./sort_bam_by_name.sh config_rnaseq.sh
```

The script will output a sorted BAM file, which can be used for generating the
count matrix using `featureCounts`. Note, this file is sorted by name, unlike
the output of the `star_alignment.sh` script, which is sorted by coordinate.

### Generating count matrix

The final step in the analysis is to generate the count matrix. This can be done
using tools such as `featureCounts`, which is a program in the `subread` package
that can be used to count the number of reads that map to each gene in the
reference genome. In the script, `featureCounts` produce the count matrix by
only considering uniquely mapped reads to increase the confidence of the
expression estimates. However, make sure to check the quality of the alignment
to ensure that uniquely mapped reads are sufficient for the analysis (should be
approximately a minimum of 70-80% of the total reads).

We can perform this via the `featureCounts_gene.sh` script.

```bash
./featureCounts_gene.sh config_rnaseq.sh
```

The script will output a count matrix, which can be used for downstream
differential expression analysis. Here, the count matrix is generated at the
gene level.

### Visualizing the RNA-seq Tracks

Visually inspecting the RNA-seq reads on different regions of the genome
is an important way to understand the expression of the genes, especially
when integrating gene expression with other genomic data. This can be done
using tools such as the Integrative Genomics Viewer (IGV), which is a
high-performance visualization tool for interactive exploration of large,
integrated genomic datasets.

Visualizing the genome track requires a compatible file format, such as the
`bigWig` or `bedgraph` format. In addition, the reads must be normalized so that 
the tracks can be compared across different samples.

The `bam_to_coverage.sh` script can be used to convert the BAM file to the bigWig
format.

Here, we use `deepTools` to convert the BAM file to the bigWig format. Ensure
this python package is installed in your environment. In addition, before
converting the BAM file to bigWig, it is important to sort the BAM file by
coordinate (which is already done as part of STAR's alignment step). We would
also need to index the BAM file using the `samtools index` command. These steps
are integrated into the `bam_to_coverage.sh` script.

```bash
./bam_to_coverage.sh config_rnaseq.sh
```

For transcriptomic data, which are splice-aware, it is important not to use the
`--extendReads` option, as this will extend the reads across the splice junctions,
which is not biologically meaningful. For choosing the normalization method, the
`CPM` (counts per million) method is recommended, as it is robust to differences
in library size and can be used to compare between samples. If relative expression
between genes is important, the `RPKM` (reads per kilobase per million) or the
`FPKM` (fragments per kilobase per million) method can be used as they take into
account the gene length.