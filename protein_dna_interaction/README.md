---
layout: default
title: Protein-DNA Interaction
nav_order: 3
---

# Protein-DNA Interaction (ChIP-seq, CUT&Tag, CUT&RUN)

One of the main mechanism of transcriptional regulation involves the interaction
between transcriptional regulators, such as transcription factors and histones,
with DNA. The study of these interactions is crucial to understand the
regulation of gene expression. Here, we describe the processing and analysis of
raw data files for ChIP-seq, CUT&Tag, and CUT&RUN experiments.

**Table of contents**

- [Protein-DNA Interaction (ChIP-seq, CUT\&Tag, CUT\&RUN)](#protein-dna-interaction-chip-seq-cuttag-cutrun)
  - [Configuration](#configuration)
  - [Pre-alignment processing](#pre-alignment-processing)
    - [Quality control](#quality-control)
    - [Trimming](#trimming)
    - [Quality control after trimming](#quality-control-after-trimming)
  - [Alignment](#alignment)
    - [Building the reference genome index](#building-the-reference-genome-index)
    - [Alignment of reads](#alignment-of-reads)
    - [Spike-in alignment (optional)](#spike-in-alignment-optional)
  - [Post-alignment QC and processing](#post-alignment-qc-and-processing)
    - [Alignment rate](#alignment-rate)
    - [Duplicate removal](#duplicate-removal)
    - [Fragment size distribution](#fragment-size-distribution)


## Configuration

Ensure that all the scripts in this folder are executable. If not, run the
following command:

```bash
chmod -R a+x ./protein_dna_interaction
```

Make sure to set up the proper parameters in the `config_chipseq.sh` file. The
[configuration file](config_chipseq.sh) is an example bash script that sets up
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
./fastqc_pre-trimming.sh config_chipseq.sh
```

### Trimming

After checking the quality of the raw data, it is important to remove low quality
bases and adapter sequences. This can be done using tools such as `fastp`.

```bash
./trimming_reads.sh config_chipseq.sh
```

### Quality control after trimming

After trimming the reads, we can again assess the quality of the trimmed reads
and examine whether the trimming procedure has improved the quality of the reads.
We should expect low quality bases are removed, along with adapter sequences.

```bash
./fastqc_post-trimming.sh config_chipseq.sh
```

<br><br/>

## Alignment

The heaviest computational step involves the alignment of the reads to the
reference genome. A common program used to align reads is `bowtie2`, which is
a fast and memory-efficient tool for aligning sequencing reads to long reference
sequences. Unlike the alignment of transcriptomic data, which is splice-aware,
reads generated from genomic DNA (either IP of fragmented chromatin or
tagmented DNA) are aligned against the entire mappable region of the genome.
Although without the need for understanding the splicing pattern, the regions of
alignment becomes more complex due to different genomic features associated with
the DNA-protein interaction, for example, the presence of repetitive elements,
heterochromatin, and methylation patterns.

The processing of epigenomic data usually requires a thorough annotation of the
genomic features of the reference genome, including the "blacklisted" regions
where the alignment is not reliable, as well as additional annotations such as
promoters and enhancers that are critical to the interpretation of the data.
Therefore, as of `2024`, the most common reference genome used for human is
`hg38` and for mouse is `mm10`. It is recommended to use them due to their
comprehensive annotations. Although the mouse genome has been updated to `mm39`,
the annotations are not as comprehensive as `mm10`.

### Building the reference genome index

Before aligning the reads, we need to build the index of the reference genome.
This is a one-time step and the index can be reused for multiple samples. The
index is a set of files that `bowtie2` uses to quickly align the reads to the
reference genome. The index is built using the `bowtie2-build` command, which
requires the reference genome in `fasta` format.

The following command builds the index of the reference genome:

```bash
bowtie2-build --threads <threads> <reference_genome.fa> <index_name>
```

Note that `<index_name>` is the prefix of the index files, rather than the name
of a directory. To keep the index files organized, it is recommended to create
a separate directory for only the bowtie2 index files. For example, the resulting
file can be stored at `/path/to/reference_genome/bowtie2_index/<index_name>`.
Index name is important for `bowtie2` to properly identify the index files when
aligning the reads.

As an example for the `mm10` mouse genome, the following command builds the index:

```bash
# Assuming mm10.fa is the fasta file of the reference genome
cd /path/to/reference_genome
mkdir ./mm10/bowtie2_index
bowtie2-build --threads  ./mm10/mm10.fa ./reference_genome/bowtie2_index/mm10
```

This will lead to a set of files in the `./mm10/bowtie2_index` directory with
the prefix `mm10` that are used by `bowtie2` for alignment. These files include
`mm10.1.bt2`, `mm10.2.bt2`, `mm10.3.bt2`, `mm10.4.bt2`, `mm10.rev.1.bt2`, and
`mm10.rev.2.bt2`.

### Alignment of reads

Upon having the index of the reference genome, we can now align the reads to the
genome. The alignment is performed using the `bowtie2` command. The following
command aligns the reads to the reference genome as part of the pipeline:

```bash
./bowtie2_align.sh config_chipseq.sh
```

Note that in the configuration file, it is necessary to set the path to the
reference genome index, which is the `bowtieIndex` parameter. This path must 
follow the following format: `/path/to/index_directory/index_name`.

In addition, the choice of the maximum fragment length is important for the
alignment. This value could be different depending on the specific experiment.
For CUT&Tag, in which nucleosomes are tagmented by the transposase, the maximum
fragment length is usually set to `700 bp`. For ChIP-seq, the maximum fragment
size distribution could be wider, and the maximum fragment length could be set
to `1200 bp`. However, consult your Bioanalyzer QC result for the library and
the method for library preparation to determine this value.

Here, we use the following alignment parameters:

```bash
--local --very-sensitive --no-mixed --no-discordant --phred33
```

Here is the explanation of the parameters used:

* `--local`: local alignment is performed
* `--very-sensitive`: the most sensitive alignment is performed (which could be computationally expensive)
* `--no-mixed`: do not match reads (paired reads) that have both unaligned and aligned parts
* `--no-discordant`: reads that do not aligns with the expected relative mate orientation or with the expected range of distances are not reported
* `--phred33`: the quality score is in phred33 format

In case users are performing sequencing with `read length = 25 bp`, there is no
need to trim the reads since adapter contents are not sequenced for `inserts > 25 bp`.
In this case, the `--local` parameter is not necessary and the alignment can be
performed with the following parameters:

```bash
--end-to-end --very-sensitive --no-mixed --no-discordant --phred33
```

During alignment, we have also converted SAM files to BAM files, which are the
binary version of SAM files. This means that BAM files are smaller and faster
to process.

### Spike-in alignment (optional)

In some experiments, spike-in controls are used to normalize the data. These
spike-in could result from additional spike-in of DNA fragments from another
species. Alternatively, in CUT&Tag experiment, bacterial DNA is carried over
during the manufacturing of the Tn5 transposase. However, our experience shows
that a very low percentage of reads are aligned to the spike-in genome, and
it may not be an optimal way to perform spike-in normalization.

To align the spike-in reads, the following command can be used:

```bash
bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant  --phred33 \
    -I 10 -X 700 \
    -p ${threads} -x ${bowtieIndex} \   # spike-in genome index
    -1 sample_R1_trimmed.fastq.gz \
    -2 sample_R2_trimmed.fastq.gz \
    -S output_bowtie2.sam &> bowtie2_summary.txt
```

<br><br/>

## Post-alignment QC and processing

After aligning the reads to the reference genome, we perform post-alignment QC
and processing.

### Alignment rate

The alignment rate is show in the `bowtie2.txt` file under the `alignment/bowtie2_summary` folder. The alignment rate
is calculated as the number of reads that are aligned to the reference genome 
(both uniquely aligned and multi-mapping reads) divided by the total number of
reads in the fastq files. It is expected that a successful alignment rate is about 80-90%.

### Duplicate removal

After alignment, we remove duplicates using `picard` tools. The removal of
duplicates is important to ensure that only unique reads are used for the
downstream analysis. This is because PCR duplicates are not informative and
will cause bias when algorithms tries to identify the enriched regions.

While many downstream analysis tools can remove duplicates on-the-fly, it is
still recommended to remove duplicates at this step to ensure that the
downstream analysis is consistent and reproducible. In this case, when
specifying some downstream analysis parameters, we need to ensure that these
downstream tools `ignore duplicates`.

The following command removes duplicates from the aligned BAM file:

```bash
./dedup_bam.sh config_chipseq.sh
```

Note that we also use `samtools` to sort the BAM file before removing duplicates.
This is required by `picard`. Furthermore, the path to the `picard.jar` file
must be set in the configuration file.

### Fragment size distribution

The fragment size distribution is important to assess the quality of the library
and to determine the size of the DNA fragments that are sequenced. This result
is very likely different from the fragment size distribution of the original
library QC by Bioanalyzer, since both the alignment step and duplicate removal
step could affect the fragment size distribution.

The SAM/BAM format stores the size of each fragment in the 9th column of the
file. Here, we use the `fragment_size.sh` script to extract the fragment size
calculate the distribution.

```bash
./fragment_size.sh config_atacseq.sh
```

The size distribution can be plotted in R, using the `ggplot2` and the `ggpubr`
packages. We have provided an R script that can be run directly from the command
line to generate the fragment size distribution plot.

To plot the fragment size distribution, run the following command:

```bash
Rscript plot_frag_size.R <frag_size_file_dir> <sample_matrix> [<plot_width>] [<plot_height>]
```

Note that although the script is similar to the one that is used for the
ATAC-seq pipeline, the usage of this script is different. By specifying the
required parameters, the script will generate multiple output files, one for
each mark used in the experiment based on the sample matrix provided. To
generate the plot, ensure that `<frag_size_file_dir>` are set to the same path
as the output of the `fragment_size.sh` script.

Furthermore, to be able to plot, a `<sample_matrix>` must have been set up.
It is up to the user to fill in the metadata of the samples in this `csv` file,
based on how the experiment was performed. To generate a baseline csv file
for all samples, users can make use of the `utils/generate_sample_matrix.py` 
script, which will create a sample matrix with the required columns. In this
case, the `plot_frag_size.R` scripts requires that the `sample_matrix` file
contains the following columns:

- `Label`: the label of the sample, which must match the file name of the
  fragment size distribution file without the `_fragment_size.txt` suffix.
- `Mark`: the protein target used. A plot will be generated for each unique
  mark in this column.

