# This is an example configuration file for the analysis of transcriptomic data
# using the RNA-seq pipeline. The configuration file is in bash format.
# The configuration file is divided into sections, each section corresponding to
# a step in the pipeline. Each section contains a list of parameters that are
# used by the pipeline to perform the analysis. The parameters are defined as
# key-value pairs, where the key is the name of the parameter and the value is
# the value of the parameter. The parameters are used to specify the input files,
# the output files, the parameters of the analysis, and the parameters of the
# software packages used in the analysis.

# For most parameters, default values are provided. The default values are used
# to establish a predefined directory structure.

# Note that the input file name must follow illumina's naming convention:
# <sample_name>_<lane>_R<read_number>_001.fastq.gz

# ====== Global parameters =====================================================
# Path to the project home directory
projPath="/your/project/path"

# Path to the directory containing the FASTQ files
fastqDir="${projPath}/fastq" 

# Lane identifier
lane="L001" 

# Number of threads to use for the analysis (check your CPU for this value)
threads="8" 

# ======= Data preprocessing ===================================================
# Path to the ORA reference file directory, required by Illumina's DRAGEN ORA
# ORA reference files can be found in https://support.illumina.com/downloads/ora-decompression-reference-files.html
# The directory should end with the species name of the reference named by Illumina, e.g. /path/to/ora_ref/mus_musculus
oraRefDir="/path/to/ora_ref"

# Path to the directory containing the FASTQC output
fastqcOutDir="${projPath}/fastqc_output"

# Path to the directory containing the trimmed FASTQ files
trimmedDir="${projPath}/fastq_trimmed" 

# Path to the Trimmomatic JAR file [Deprecated]
# Ignored if use fastp for trimming
TRIMMOMATIC_JAR="/path/to/trimmomatic.jar" 

# Path to the Trimmomatic adapter file [Deprecated]
# Ignored if use fastp for trimming
TRIMMOMATIC_ADAPTER_FILE="/path/to/adapter.fa" 
# Note that the adapter must be chosen corresponding to the library generation method
# These fasta files are usually stored along with the Trimmomatic JAR file

### For checking rRNA contamination

# Whether rRNA filtering is enabled
# If set to "true", the pipeline will filter out reads that align to rRNA
# If set to "false", the pipeline will skip the rRNA filtering step
rRNAFiltering="false" # Must be either "true" or "false"

# Path to the rRNA FASTA file
# This is a species specific file containing the rRNA sequences
rRNAFastaFile="/path/to/rRNA.fa"

# Path to the rRNA GTF file
rRNAGtfFile="/path/to/rRNA.gtf"

# Path to the SortMeRNA reference file
# SortMeRNA already provides some pre-built databases that include multi-taxon rRNA species
sortmernaRef="/path/to/sortmerna_ref"

# SortmeRNA working directory
sortmernaWorkDir="${USER}/sortmerna"

# Path to the output directory for the rRNA alignment, fastq files that align to rRNA
rRNAalignDir="${projPath}/rRNA_filtering/rRNA_aligned"

# Path to the output directory for the rRNA contamination free fastq files
rRNAcleanDir="${projPath}/rRNA_filtering/rRNA_cleaned"

# ====== Alignment =============================================================
# Path to the STAR index directory
STARIndexDir="/path/to/STARIndex"
# if you do not have a STAR index, the pipeline will create one in this directory

# Path to the genome FASTA file
genomeFastaFile="/path/to/genome.fa"

# Path to the genome GTF file
genomeGtfFile="/path/to/genome.gtf" 
# Note that the genome FASTA and GTF files are not used if you already have a STAR index

# Path to the output directory for the alignment
STAROutputDir="${projPath}/star_output"

# Temporary FIFO file directory
# STAR requires a temporary directory to store FIFO files, which should be on
# a linux filesystem.
# If useTempFIFO is set to "false", this temporary directory will the same
# as the STAROutputDir
useTempFIFO="false"   # Must be either "true" or "false"
tempFIFODir="/tmp"

### Only used when running star_alignment_short_fragment.sh
# Reduce requirement for filtering based on fragment length
# The following parameters are used to filter out reads based on fragment length:
# outFilterScoreMinOverLread: minimum score of a read alignment as a fraction of the read length
# outFilterMatchNminOverLread: minimum number of matched bases in a read alignment as a fraction of the read length
# outFilterMismatchNmax: maximum number of mismatches in a read alignment
filterLreadScore="0"
filterLreadMatch="0"
filterMismatchNmax="2"

# ====== Post-alignment processing =============================================
# Path to the directory containing the BAM files sorted by name
nameSortedDir="${projPath}/samtools_sort_name"

# Path to the output directory of featureCounts
featureCountsDir="${projPath}/feature_counts"

# Name of the output file for the featureCounts table
outCountFile="raw_counts.txt"

# Path to the coverage files for visualizing RNAseq tracks
coverageFilePath="${projPath}/coverage_tracks"

# Coverage file format, either "bedgraph" or "bigwig"
outCoverageFormat="bigwig"

# Normalization method for the coverage files, one of RPKM, CPM, BPM, RPGC, None
coverageNorm="CPM"

# Effective genome size for the coverage files, used for normalization
# Refer to utils/resources
genomeSize="XXXXXXXXXX"

# Bin size for the coverage files
binSize="20"




# [END] 