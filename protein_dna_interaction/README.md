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
    - [Standardizing the sample names](#standardizing-the-sample-names)
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
    - [MACS2 peak calling](#macs2-peak-calling)
      - [Annotating samples for peak calling](#annotating-samples-for-peak-calling)
      - [Extracting samples for peak calling type](#extracting-samples-for-peak-calling-type)
      - [Peak calling](#peak-calling-1)
    - [goPeak](#gopeak)
  - [Processing called peaks](#processing-called-peaks)
    - [Blacklist filtering](#blacklist-filtering)
    - [Peak statistics](#peak-statistics)
    - [Visualizing coverage tracks](#visualizing-coverage-tracks)
    - [Peak coverage and distribution](#peak-coverage-and-distribution)
  - [Identifying consensus peaks](#identifying-consensus-peaks)
    - [Calculating replicate coverage over peaks](#calculating-replicate-coverage-over-peaks)
    - [Intersecting peaks](#intersecting-peaks)


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

### Standardizing the sample names

Immediately after alignment, it is recommended to standardize the sample names
to ensure that the downstream analysis clear. The sample names, as well as the
file names, should clearly indicate the unique sample identifier, the mark
used, and the replicate number.

Regardless, a sample matrix csv file is required for the pipeline to run
effectively, since the downstream analysis scripts will process samples by
grouping them based on the mark used.

To initialize a sample matrix csv, users can use the `utils/generate_sample_matrix.py`. See the [documentation](https://j-y26.github.io/Epigenomic-Pipelines/utils/docs/generate_sample_matrix.html) for more information. This script will generate a baseline sample matrix with the required columns. Users can then fill in the metadata of the samples in this `csv` file, based on how the experiment was performed.

The csv file must contain the following columns:

- `Label`: the label of the sample, which must match the file name of the fastq files
- `Mark`: the protein target used

Users should manually fill in the metadata of the samples in this `csv` file,
based on how the experiment was performed. The `Label` column should match the
desired sample name, and the "Sample" column will match the original file name
used earlier in the pipeline.

The pipeline has also provided an easy way to standardize the file names
based on the sample matrix. The python script `utils/file_name_conversion.py`
can be used to convert the file names based on the sample matrix, from the
`Sample` column to the `Label` column. See the [documentation](https://j-y26.github.io/Epigenomic-Pipelines/utils/docs/file_name_conversion.html) for more information.

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
./filter_index_bam.sh config_chipseq.sh
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
./fragment_size.sh config_chipseq.sh
```

The size distribution can be plotted in R, using the `ggplot2` and the `ggpubr`
packages. We have provided an R script that can be run directly from the command
line to generate the fragment size distribution plot.

To plot the fragment size distribution, run the following command:

```bash
Rscript plot_frag_size.R <frag_size_file_dir> <sample_matrix> <max_fragment_length> [<plot_width>] [<plot_height>]
```

Note that although the script is similar to the one that is used for the
ATAC-seq pipeline, the usage of this script is different. By specifying the
required parameters, the script will generate multiple output files, one for
each mark used in the experiment based on the sample matrix provided. To
generate the plot, ensure that `<frag_size_file_dir>` are set to the same path
as the output of the `fragment_size.sh` script.

Furthermore, to be able to plot, a `<sample_matrix>` must have been set up.
It is up to the user to fill in the metadata of the samples in this `csv` file,
based on how the experiment was performed. The `plot_frag_size.R` scripts
requires that the `sample_matrix` file contains the following columns:

- `Label`: the label of the sample, which must match the file name of the
  fragment size distribution file without the `_fragment_size.txt` suffix.
- `Mark`: the protein target used. A plot will be generated for each unique
  mark in this column.

The `<max_fragment_length>` parameter is used to specify the maximum fragment
length to be plotted. This value is used to set the x-axis limit of the plot.
It is recommended to set this value to the same as the `maxInsertLength` value
set in the `config_chipseq.sh` file to capture all the fragments. This is a
required input parameter.

The `<plot_width>` and `<plot_height>` parameters are optional and are used to
specify the width and height of the plot. The default values are `8` and `6`,
respectively.

An output PDF figure will be saved in the same directory as the input file.

An example of a streamlined analysis after running the `fragment_size.sh` script
is shown below:

```bash
Rscript plot_frag_size.R ${alignmentDir}/fragment_size sample_matrix.csv ${maxInsertLength}
```

#### Fragment size QC by deepTools

Alternatively, `deepTools` can be used to generate the fragment size distribution
plot and to assess the quality of the library. Note that to use the following
script in the pipeline, only paired-end data is supported.

Before running the script, we need to first extract the samples corresponding to
each mark from the sample matrix. To this end, we use the python script
`samples_for_mark.py`:

```bash
python3 samples_for_mark.py <sample_matrix.csv> <output_dir>
```

Here, `<sample_matrix.csv>` is the sample matrix file, which must contains a
`Label` and a `Mark` column. If you have not set up the sample matrix, you can
use the `utils/generate_sample_matrix.py` script to generate a baseline sample
matrix as described [here](https://j-y26.github.io/Epigenomic-Pipelines/utils/docs/generate_sample_matrix.html).

An working example here is:

```bash
python3 samples_for_mark.py sample_matrix.csv ${markedSamples}
```

After extracting the samples, we can run the `fragment_size_deeptools.sh` script
to generate the fragment size distribution plot using `deepTools`. The script
requires the following parameters:

- `markedSamples`: the path to the directory containing the samples for each mark (as generated by the `samples_for_mark.py` script)
- `bamSuffix`: the suffix of bam files where by adding the suffix to the sample name, the file name of the bam file can be obtained. For example, for a sample named `sample1`, the bam file name is `sample1_nodup.bam` if the suffix is set to `_nodup.bam`.

The following command runs the script:

```bash
./fragment_size_deeptools.sh config_chipseq.sh
```

The script will generate a fragment size distribution plot for each mark in the
 `fragment_size` folder
under the `alignment` directory. The plot for each mark is saved as a PDF file, which is a
histogram of the fragment size distribution. In addition, a text file is also
generated that contains the fragment size distribution data.

Unlike the manual extraction of the fragment size from the SAM/BAM files, we
take advantage of `deepTools's` feature for removing the blacklisted regions.
Hence, a `BED` file containing the blacklisted regions is required. The path to
the blacklisted regions file must be set in the configuration file. For common
genomes, our pipeline has already included the blacklisted regions file, which
is located in the [`utils` directory](https://j-y26.github.io/Epigenomic-Pipelines/utils/docs/resources.html).

### Reproducibility and coverage

The reproducibility of the ChIP-seq data is assessed by the correlation of the
read coverage between replicates. The read coverage is calculated using the
`deepTools` package. We will generate a coverage plot and a correlation plot
to assess the quality of the mapped reads.

#### Coverage

We can first assess the coverage of the aligned reads by inspecting a coverage
plot. We can plot using the `plotCoverage` command in `deepTools`. The following
command generates the coverage plot:

```bash
./coveragePlot.sh config_chipseq.sh
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

Both plots requires computing the coverage of the aligned reads. Since we are
interested in the quality of mapped reads for the entire genome, we will use
the `bins` mode in the `multiBamSummary` command to calculate the coverage of
the entire genome. The `bins` mode calculates the coverage of the entire genome
in non-overlapping bins. The default bin size is 10,000 bp, which is
recommended for the entire genome. Note that users can specify the bin size
in the configuration file. However, increasing the resolution (smaller bin size)
will increase the computational time and memory usage, and vice versa.

Note, if you have not yet set up the sample matrix and extracting the samples
for each mark, you can follow the description in the section 
[Fragment size QC by deepTools](#fragment-size-qc-by-deeptools) to set up the
sample matrix and extract the samples for each mark.

The following command calculates the coverage of the aligned reads:

```bash
./multiBamSummary.sh config_chipseq.sh
```

This will generate a computed coverage file for all bam files assigned to each
mark in the `bam_qc` folder under the `alignment` directory.

#### Reproducibility

The reproducibility of the ChIP-seq data is assessed by the correlation of the
read coverage between replicates. The correlation plot is generated using the
`plotCorrelation` command in `deepTools`. Furthermore, we can assess the
reproducibility by performing a principal component analysis (PCA) of the
read coverage. The PCA plot is generated using the `plotPCA` command in
`deepTools`.

Using the following script, we together generate the correlation and PCA plots:

```bash
./bamQCPlots.sh config_chipseq.sh
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
choose to increase or decrease this value. When mapping regions of a genome
that is targeted by a specific protein, it is expected that most parts of the
genome are not covered by the reads. In this case, the `ntop` parameter can be
set to a lower value to capture the variance of the regions that are covered
by the reads. Yet, a proper value depends on the specific mark used and the
quality of the antibody used in the experiment.

### Mitochondrial reads

The presence of mitochondrial reads is a common issue in genomic data analysis.
High levels of mitochondrial reads can be indicative of poor quality of the
library, and can be a result of contamination during library preparation or
sequencing. Therefore, it is important to assess the percentage of mitochondrial
reads in the data.

The following command calculate the mitochondrial and total reads counts
and generates a csv file.

```bash
./mt_reads.sh config_chipseq.sh
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

After the post-alignment QC and processing, we can now identify the enriched
regions in the genome. ChIP-seq, CUT&Tag, and CUT&RUN experiments are designed
to obtain reads localized around protein binding sites. Regions that are truly
bound by the protein of interest are expected to have a higher read coverage
compared to the background, termed enriched regions or peaks. The identification
of these peaks is performed using peak calling algorithms.

### MACS2 peak calling

Here, we use the `MACS2` peak caller to identify the enriched regions in the
genome. `MACS2` is a popular peak caller that is widely used in the field of
genomic data analysis. It is specifically designed for identifying enriched
regions in ChIP-seq data.

It is important to note that the choice of the parameters in `MACS2` is critical to obtain correctly called peaks. Specifically, CUT&Tag and ChIP-seq data have different characteristics, and the parameters used for peak calling should be adjusted accordingly. Furthermore, depending on the distribution
characteristics of the mark, narrow or broad peaks should be called independently.

#### Annotating samples for peak calling

Refer to the documentation on histone/chromatin mark peak profiles for more information on the distribution profiles of some common histone marks and RNA polymerase marks. A summary of the distribution profiles of some common histone marks and RNA polymerase marks is provided below:

[Summary of Histone/Chromatin Peak Profiles](../utils/docs/histone_peak_profile.md)

Note that in the `sample_matrix` file, a column titled `Peak_Type` is required to specify the type of peak calling to be performed. The `Peak_Type` column should contain either `narrow` or `broad` for each sample. This column is used to specify the type of peak calling to be performed by `MACS2`.

#### Extracting samples for peak calling type

For uniformed processing of peak calling, we want to extract the list of samples (files) corresponding to each type of peaks. This is done by the `extract_samples_by_peak_type.py` script. The script will be invoked as follows:

```bash
python3 extract_samples_by_peak_type.py <sample_matrix> <output_dir>
```

It is recommended to extract the samples for each peak type into a file into
the `peak_calling` directory. An example of the script is shown below:

```bash
python3 extract_samples_by_peak_type.py sample_matrix.csv ${peakCallingDir}
```

#### Peak calling

Thorough consideration of the parameters used in peak calling is
required. The parameters used in the `macs2 callpeak` command are critical to
accurately identify the enriched open chromatin regions. The following parameters
are used in the `macs2 callpeak` command:

- `-f BAMPE`: the format of the input file
- `-g ${genomeSize}`: the effective genome size for the reference genome. This number can be determined from deepTools. Some common effective genome sizes are listed along with this [pipeline](../utils/docs/resources.md).
- `-q ${qValue}`: the q-value threshold for peak calling. This is the minimum FDR at which a peak is called significant. The default value is 0.05, but users can choose to increase or decrease this value. As a common practice for ChIP-seq, a q-value threshold of 0.01 is recommended as a starting point.
- `-p ${pValue}`: the p-value threshold for peak calling. Unlike the `q-value`, the `p-value` is not corrected for multiple testing and represents the probability of observing a peak by chance. If the `p-value` is set, the `q-value` is ignored. We do not recommend using the `p-value` for peak calling since it increases the false positive rate.
- `-B`: generate bedGraph files of the pileup, which can be used to visualize the pileup in a genome browser.
- `--SPMR`: use the signal per million reads (SPMR) as the normalization method. This does not affect the peak calling, but is used to normalize the signal for visualization.
- `--call-summits`: call the summits of the peaks. This is recommended for ChIP-seq data, as it identifies the exact location of the peak.

As mentioned earlier, a cutoff is required to filter the peaks to keep only the
significant ones. The `q-value` is the most commonly used cutoff, and it is
recommended to set the `q-value` to `0.01` as a starting point. However, the
optimal cutoff value might not be the same for all experiments, and it is
may be necessary to run cutoff analysis to determine the optimal cutoff value.
`MACS2` provides a cuttoff analysis mode to determine the optimal cutoff value.
To do so, in the configuration file, set the `cutoffAnalysis` parameter to `true`.

More specifically, users should consider the following parameters:

- `--nomodel`: whether to disable the shifting model.
- `--extsize`: the extension size of the reads. This is the size of the DNA fragment that is sequenced. The extension size should be set to the average fragment size of the library. It is required is no shifting model is used.
- `--nolambda`: whether to disable the local lambda calculation.

The following command performs peak calling using `MACS2`:

```bash
./macs2_peak_calling.sh config_chipseq.sh
```

### goPeak

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
./blacklist_filtering.sh config_chipseq.sh
```

### Peak statistics

The statistics of the called peaks are calculated to assess the quality of the
peaks and provide a summary of the peak calling results. The statistics include
the number of peaks called and the distribution of the peak widths.

The following command calculates the peak statistics:

```bash
Rscript plot_peak_stats.R <peak_bed_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]
```

The `<peak_bed_dir>` is the directory containing the filtered peak bed files.
The `<output_dir>` is the directory where the peak statistics plots and data
will be saved. The `<plot_width>` and `<plot_height>` parameters are optional
and are used to specify the width and height of the plot. The default values
are `8` and `6`, respectively.

An example of a streamlined analysis after running the `blacklist_filtering.sh`
script is shown below:

```bash
peakStatsDir=${peakCallingDir}/peak_stats/${rawPeaks}
filteredPeaks=${peakCallingDir}/filtered_peaks/${rawPeaks}
mkdir -p ${peakStatsDir}
Rscript plot_peak_stats.R ${filteredPeaks} sample_matrix.csv ${peakStatsDir}
```

Note that while most peaks usually have a width less than 1 kb, some peaks may
be much wider. Running the script generates a summary statistics that allows users
to inspect the peaks size distribution in its naive form. On the other hand,
to allow better comparison between samples, if the maximum peak width of a sample
is greater than 10 kb, the script will only consider the maximum peak width of
at the 90th percentile across all samples of the same mark.
This is to ensure that the peak width
distribution is not skewed by a few samples with very wide peaks.

### Visualizing coverage tracks

The coverage of the aligned reads can be visualized in a genome browser to
inspect the quality of the data and the called peaks. We will convert the `BAM`
files used for peak calling to `bigWig` or `bedgraph` files, which are used to 
visualize the coverage in a genome browser. The files are generated using the
`bamCoverage` command in `deepTools`.

The following command generates the `bigWig` or `bedgraph` files:

```bash
./bam_to_coverage.sh config_chipseq.sh
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
./computeMatrix.sh config_chipseq.sh
```

To generate the profile and heatmap of the coverage, run the following command:

```bash
./plotProfileHeatmap.sh config_chipseq.sh
```

## Identifying consensus peaks

To allow confident peak detection, it is essential to identify the consensus
peaks across replicates. The consensus peaks are the peaks that are called in
multiple replicates, and are considered to be the true peaks.

While there are many exiting methods that can be used to identify the consensus
peaks, they usually involve identifying regions of peaks that are called in
multiple replicates. Here, we provide two ways to identify the consensus peaks.

Before calling any peaks, it is important to identify the samples that are
considered as replicates. Here, we have provided a python script that can be
used to identify the replicates based on the sample matrix.

The following command identifies the replicates based on the sample matrix:

```bash
python3 extract_replicated_peaks.py <sample_matrix> <output_dir>
```

In the following scripts we require that these replicate information must be
placed directly under the `peak_calling/consensus_peaks` directory. Therefore,
run the following script to ensure the replicate identification files are placed
in the correct directory:

```bash
consensusPeaksDir=${peakCallingDir}/consensus_peaks
python3 extract_replicated_peaks.py sample_matrix.csv ${consensusPeaksDir}
```

### Calculating replicate coverage over peaks

In this method, we first categorize the genome in to non-equal0size bins based on
the peak overlap between different replicates. The coverage (number of replicates
with peaks) covering each bin is then calculated to generate a bedGraph file.
This bedGraph file can be utilized by `macs2 bdgpeakcall` to call the consensus
peaks.

The following command calculates the replicate coverage over peaks:

```bash
./find_consensus_peaks_coverage.sh config_chipseq.sh
```

Note that just as the earlier steps, identifying the correct peaks depends on
setting the correct `${rawPeaks}` variable in the configuration file. This
variable should be able to identify the directory name specifying the original
peak calling configuration.

Furthermore, an important information is the size of each chomosome. A
`genome.chrom.sizes` file is required to run the `macs2 bdgpeakcall` command.
It is a tab-separated file containing the name of the chromosome and the size
of the chromosome. Note that the chromosome names in the `genome.chrom.sizes`
file should match the chromosome names in the peak bed files and covers all
chromosomes listed in the peak bed files. See the 
[documentation](../utils/docs/update_contig_names.md) for more information.
If any issue occurs with chromosome name inconsistency.

### Intersecting peaks

Another methods is to intersect the peaks called in different replicates to
identify the consensus peaks. The `bedtools intersect` command can be used to
intersect the peaks called in different replicates.

The following command intersects the peaks called in different replicates:

```bash
./find_consensus_peaks_overlap.sh config_chipseq.sh
```

Note that it is important to define the minimum fraction of overlap between
the peaks to be considered as the same peak. This value can be set in the
configuration file. The default value is `0.5`, which means that at least 50%
of the peak must overlap to be considered as the same peak between replicates.