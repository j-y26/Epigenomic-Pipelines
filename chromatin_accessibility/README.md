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
  - [Alignment](#alignment)
    - [Building the reference genome index](#building-the-reference-genome-index)
    - [Alignment of reads](#alignment-of-reads)
    - [Spike-in alignment (optional)](#spike-in-alignment-optional)
  - [Post-alignment QC and processing](#post-alignment-qc-and-processing)
    - [Alignment rate](#alignment-rate)
    - [Duplicate removal](#duplicate-removal)
    - [Read filtering and indexing](#read-filtering-and-indexing)
    - [Fragment size distribution](#fragment-size-distribution)
      - [Extracting fragment size from SAM/BAM files](#extracting-fragment-size-from-sambam-files)
      - [Fragment size QC by deepTools](#fragment-size-qc-by-deeptools)
    - [Reproducibility and coverage](#reproducibility-and-coverage)
      - [Coverage](#coverage)
      - [Genome-wide read coverage](#genome-wide-read-coverage)
      - [Reproducibility](#reproducibility)
    - [Mitochondrial reads](#mitochondrial-reads)
  - [Peak calling](#peak-calling)
  - [Processing called peaks](#processing-called-peaks)
    - [Blacklist filtering](#blacklist-filtering)
    - [Visualizing coverage tracks](#visualizing-coverage-tracks)
    - [Peak coverage and distribution](#peak-coverage-and-distribution)

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

<br><br/>

## Alignment

The heaviest computational step involves the alignment of the reads to the
reference genome. A common program used to align reads is `bowtie2`, which is
a fast and memory-efficient tool for aligning sequencing reads to long reference
sequences. Unlike the alignment of transcriptomic data, which is splice-aware,
reads generated from tagmented genomic DNA by the Tn5 transposase are aligned 
against the entire mappable region of the genome.
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
./bowtie2_align.sh config_atacseq.sh
```

Note that in the configuration file, it is necessary to set the path to the
reference genome index, which is the `bowtieIndex` parameter. This path must 
follow the following format: `/path/to/index_directory/index_name`.

In addition, the choice of the maximum fragment length is important for the
alignment. This value could be different depending on the specific experiment.
The transposase accessible regions are usually around `200-300 bp` in length, but
depending on the tagmentation efficiency, the fragment length could be longer.
For ATAC-seq, the maximum fragment length is usually set to `1000 bp`.
However, consult your Bioanalyzer QC result for the library and
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
./dedup_bam.sh config_atacseq.sh
```

Note that we also use `samtools` to sort the BAM file before removing duplicates.
This is required by `picard`. Furthermore, the path to the `picard.jar` file
must be set in the configuration file.

### Read filtering and indexing

After removing duplicates, we filter the reads to remove read pairs that are
not properly mapped. Specifically, we use the `SAM FLAG` to filter the reads
based on the following criteria, using the `samtools` command:

- Filter out unmapped reads: `samtools view -F 0x04`
- Keep only properly paired reads: `samtools view -f 0x02`

Further indexing of the BAM file is also performed using `samtools index`.
Indexing the reads is required for downstream analysis, such as quality control
and peak calling.

The following command filters the reads and indexes the BAM file:

```bash
./filter_index_bam.sh config_atacseq.sh
```

### Fragment size distribution

#### Extracting fragment size from SAM/BAM files

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
Rscript plot_frag_size.R <frag_size_file_dir> <max_fragment_length> [<plot_width>] [<plot_height>]
```

Note that although the script is similar to the one that is used for the
ChIP-seq pipeline, the usage of this script is different. By specifying the
required parameters, the script will generate a single output file that contains
the fragment size distribution plot for all samples in the given directory. To
generate the plot, ensure that `<frag_size_file_dir>` are set to the same path
as the output of the `fragment_size.sh` script.

The `<max_fragment_length>` parameter is used to specify the maximum fragment
length to be plotted. This value is used to set the x-axis limit of the plot.
It is recommended to set this value to the same as the `maxInsertLength` value
set in the `config_atacseq.sh` file to capture all the fragments. This is a
required input parameter.

The `<plot_width>` and `<plot_height>` parameters are optional and are used to
specify the width and height of the plot. The default values are `8` and `6`,
respectively.

An output PDF figure will be saved in the same directory as the input file.

An example of a streamlined analysis after running the `fragment_size.sh` script
is shown below:

```bash
Rscript plot_frag_size.R ${alignmentDir}/fragment_size ${maxInsertLength}
```

#### Fragment size QC by deepTools

Alternatively, `deepTools` can be used to generate the fragment size distribution
plot and to assess the quality of the library. Note that to use the following
script in the pipeline, only paired-end data is supported.

The following command generates the fragment size distribution plot using
`deepTools`:

```bash
./fragment_size_deeptools.sh config_atacseq.sh
```

The fragment size distribution plot is generated in the `fragment_size` folder
under the `alignment` directory. The plot is saved as a PDF file, which is a
histogram of the fragment size distribution. In addition, a text file is also
generated that contains the fragment size distribution data.

Unlike the manual extraction of the fragment size from the SAM/BAM files, we
take advantage of `deepTools's` feature for removing the blacklisted regions.
Hence, a `BED` file containing the blacklisted regions is required. The path to
the blacklisted regions file must be set in the configuration file. For common
genomes, our pipeline has already included the blacklisted regions file, which
is located in the [`utils` directory](https://j-y26.github.io/Epigenomic-Pipelines/utils/docs/resources.html).

### Reproducibility and coverage

The reproducibility of the ATAC-seq data is assessed by the correlation of the
read coverage between replicates. The read coverage is calculated using the
`deepTools` package. We will generate a coverage plot, a correlation plot,
and a PCA plot
to assess the quality of the mapped reads.

#### Coverage

We can first assess the coverage of the aligned reads by inspecting a coverage
plot. We can plot using the `plotCoverage` command in `deepTools`. The following
command generates the coverage plot:

```bash
./coveragePlot.sh config_atacseq.sh
```

Note that here we need to specify the number of 1-bp regions to be sampled
by the program to calculate the coverage. The default value is 1,000,000, which
is recommended for the entire genome. However, users can specify the number of
regions to be sampled in the configuration file. Increasing the number of
regions will increase the computational time and memory usage, and vice versa.
In addition, here we do not sample the blacklisted regions. Thus, users are
expected to provide a path to the blacklisted regions BED file in the
configuration file.

#### Genome-wide read coverage

Both the correlation and PCA plots require computing the coverage of the aligned reads. Since we are
interested in the quality of mapped reads for the entire genome, we will use
the `bins` mode in the `multiBamSummary` command to calculate the coverage of
the entire genome. The `bins` mode calculates the coverage of the entire genome
in non-overlapping bins. The default bin size is 10,000 bp, which is
recommended for the entire genome. Note that users can specify the bin size
in the configuration file. However, increasing the resolution (smaller bin size)
will increase the computational time and memory usage, and vice versa.

The following command calculates the coverage of the aligned reads:

```bash
./multiBamSummary.sh config_atacseq.sh
```

This will generate a computed coverage file in the `bam_qc` folder under the 
`alignment` directory.

#### Reproducibility

The reproducibility of the ATAC-seq data is assessed by the correlation of the
read coverage between replicates. The correlation plot is generated using the
`plotCorrelation` command in `deepTools`. Furthermore, we can assess the
reproducibility by performing a principal component analysis (PCA) of the
read coverage. The PCA plot is generated using the `plotPCA` command in
`deepTools`.

Using the following script, we together generate the correlation and PCA plots:

```bash
./bamQCPlots.sh config_atacseq.sh
```

Note that `plotCorrelation` and `plotPCA` require the computed coverage file
generated by the `multiBamSummary` command. This means that the `multiBamSummary.sh`
script must be run before running the `bamQCPlots.sh` script.

For the correlation plot, a method used to calculate the correlation coefficient
should be specified in the configuration file. The default method is `pearson`,
which is recommended for most cases. Another method that can be used is `spearman`.

For the PCA plot, by default, not all bins with variance are used to calculate
the PCA. The `ntop` parameter specifies the number of bins with the highest
variance to be used for the PCA. The default value is 1000, but users can
choose to increase or decrease this value. Since for ATAC-seq data, heterochromatic
regions are not expected to contain many reads, a major portion of the genome,
and thus bins, will have low variance. Therefore, it is recommended to set
`ntop` to capture the highly variable regions of the genome.

### Mitochondrial reads

The presence of mitochondrial reads is a common issue in genomic data analysis.
High levels of mitochondrial reads can be indicative of poor quality of the
library, and can be a result of contamination during library preparation or
sequencing. Therefore, it is important to assess the percentage of mitochondrial
reads in the data.

The following command calculate the mitochondrial and total reads counts
and generates a csv file.

```bash
./mt_reads.sh config_atacseq.sh
```

The csv file have the following structure:

| sample | mt_reads | total_reads |
|--------|----------|-------------|
| sample1 | xxx | xxx |
| sample2 | xxx | xxx |

Based on the given information, we can easily calculate the percentage of
mitochondrial reads for each sample.

Alternatively, we can also generate a plot to visualize the percentage of
mitochondrial reads for each sample. The following command generates the plot:

```bash
python3 plot_mt_reads.py <csv file> <output_dir>
```

The `<csv file>` is the csv file generated by the `mt_reads.sh` script. The
`<output_dir>` is the directory where the plot will be saved. It is recommended
to set the output directory to the `mt_reads` folder under the `alignment` directory.

```bash
python3 plot_mt_reads.py ${alignmentDir}/mt_reads/mt_reads.csv ${alignmentDir}/mt_reads
```

<br><br/>

## Peak calling

The identification of chromatin accessible regions is performed by identifying
the fragment pileup regions, which are indicative of open chromatin regions.
This is done by calling peaks using some of the most popular peak calling
algorithms, such as `MACS2`. The peak calling step is the most important step
in the ATAC-seq analysis, as it identifies the regions of the genome that are
enriched for open chromatin.

The following command calls peaks using `MACS2`:

```bash
./macs2_call_peaks.sh config_atacseq.sh
```

The `macs2_call_peaks.sh` script calls peaks using the `macs2 callpeak` command.

However, thorough consideration of the parameters used in peak calling is
required. The parameters used in the `macs2 callpeak` command are critical to
accurately identify the enriched open chromatin regions. The following parameters
are used in the `macs2 callpeak` command:

- `-f BAMPE`: the format of the input file
- `-g ${genomeSize}`: the effective genome size for the reference genome. This number can be determined from deepTools. Some common effective genome sizes are listed along with this [pipeline](../utils/docs/resources.md).
- `-q ${qValue}`: the q-value threshold for peak calling. This is the minimum FDR at which a peak is called significant. The default value is 0.05, but users can choose to increase or decrease this value. As a common practice for ATAC-seq, a q-value threshold of 0.01 is recommended as a starting point.
- `-p ${pValue}`: the p-value threshold for peak calling. Unlike the `q-value`, the `p-value` is not corrected for multiple testing and represents the probability of observing a peak by chance. If the `p-value` is set, the `q-value` is ignored. We do not recommend using the `p-value` for peak calling since it increases the false positive rate.
- `--nolambda`: do not calculate the local lambda, which is used to estimate the background noise. This is recommended for ATAC-seq data since no control (input) sample is used in an ATAC-seq experiment.
- `-B`: generate bedGraph files of the pileup, which can be used to visualize the pileup in a genome browser.
- `--SPMR`: use the signal per million reads (SPMR) as the normalization method. This does not affect the peak calling, but is used to normalize the signal for visualization.
- `--call-summits`: call the summits of the peaks. This is recommended for ATAC-seq data, as it identifies the exact location of the peak.

Due to the cumulative nature of ATAC-seq peaks present in accessible promoters and enhancers, a narrow peak is expected. Therefore, the `MACS2` algorithm here operates on its default setting to call narrow peaks.

As mentioned earlier, a cutoff is required to filter the peaks to keep only the
significant ones. The `q-value` is the most commonly used cutoff, and it is
recommended to set the `q-value` to `0.01` as a starting point. However, the
optimal cutoff value might not be the same for all experiments, and it is
may be necessary to run cutoff analysis to determine the optimal cutoff value.
`MACS2` provides a cuttoff analysis mode to determine the optimal cutoff value.
To do so, in the configuration file, set the `cutoffAnalysis` parameter to `true`.

## Processing called peaks

### Blacklist filtering

The called peaks are filtered to remove the blacklisted regions. The blacklisted
regions are regions of the genome where the alignment is not reliable, and are
commonly used to filter out false positive peaks. Therefore, we do not want to
analyze the peaks that are called in the blacklisted regions.

Here, we use the `bedtools` command to remove the `.narrowPeak` files output by
`MACS2` that overlap with the blacklisted regions.

`bedtools` requires that the naming convention of the chromosomes are consistent
between the `.narrowPeak` file and the blacklisted regions file. Depending on
the reference genome used, the naming convention of the chromosomes may be
different. Here, we provide an easy script to remove the leading `chr` from the
chromosome names in a `bed` file.

```bash
sed 's/^chr//g' blacklist_genebank_naming.bed > blacklist_ensembl_naming.bed
```

The following command filters the called peaks to remove the blacklisted regions:

```bash
./blacklist_filtering.sh config_atacseq.sh
```

### Visualizing coverage tracks

The coverage of the aligned reads can be visualized in a genome browser to
inspect the quality of the data and the called peaks. We will convert the `BAM`
files used for peak calling to `bigWig` or `bedgraph` files, which are used to 
visualize the coverage in a genome browser. The files are generated using the
`bamCoverage` command in `deepTools`.

The following command generates the `bigWig` or `bedgraph` files:

```bash
./bam_to_coverage.sh config_atacseq.sh
```

### Peak coverage and distribution

The coverage of the called peaks is assessed to ensure that the peaks are
enriched for open chromatin. The coverage of the peaks is calculated using the
`computeMatrix`, and the distribution of the coverage is visualized using the
`plotProfile` and `plotHeatmap` commands in `deepTools`.

To do so, we will first compute the matrix of the coverage of the peaks, and
then plot the profile and heatmap of the coverage.

The following command computes the matrix of the coverage of the peaks:

```bash
./computeMatrix.sh config_atacseq.sh
```