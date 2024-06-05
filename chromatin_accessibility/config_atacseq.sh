# This is an example configuration file for the analysis of transcriptomic data
# using the ATAC-seq pipeline. The configuration file is in bash format.
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
projPath="/your/project/path"   # <-- Update this

# Path to the directory containing the FASTQ files
fastqDir="${projPath}/fastq" 

# Lane identifier
lane="L001"   # <-- Update this

# Number of threads to use for the analysis (check your CPU for this value)
threads="8"   # <-- Update this

# ======= Data preprocessing ===================================================
# Path to the directory containing the FASTQC output
fastqcOutDir="${projPath}/fastqc_output" 

# Path to the directory containing the trimmed FASTQ files
trimmedDir="${projPath}/fastq_trimmed" 

# ====== Alignment =============================================================
# Path to the Bowtie2 index
# Note that the end of the below path must be the prefix of the index
bowtieIndex="/path/to/bowtie2_index"

# Path to the genome FASTA file
genomeFastaFile="/path/to/genome.fa"  # <-- Update this

# Path to the genome GTF file
genomeGtfFile="/path/to/genome.gtf"   # <-- Update this

# Path to the output directory for the alignment
alignmentDir="${projPath}/alignment"

# Maximum insert length for paired-end reads
maxInsertLength="1000"

# Whether to preserve the human-readable SAM files, must be their "true" or "false"
keepSAM="false"

# ====== Post-alignment processing =============================================
# Path to Picard JAR file
PICARD_JAR="/path/to/picard.jar"  # <-- Update this

# Path to the blacklist file for the specified genome
blacklistFile="/path/to/blacklist.bed"  # <-- Update this

# bin size for the multiBamSummary
binSize="10000"

####  Coverage plot
# Number of 1bp regions to sample for the coverage plot
numberOfSamples="1000000"

# Height of the coverage plot
coveragePlotHeight="5"

# Width of the coverage plot
coveragePlotWidth="15"

#### Correlation plot
# Method for calculating the correlation, either "pearson" or "spearman"
correlationMethod="pearson"

# Color map for the correlation plot
colorMap="viridis"

# Correlation plot height
correlationPlotHeight="9.5"

# Correlation plot width
correlationPlotWidth="11"

#### PCA plot
# Number of top variance bins used for PCA plot
ntop="1000"

# PCA plot height
pcaPlotHeight="10"

# PCA plot width
pcaPlotWidth="10"

#### Mitochondrial reads
# Chromosome name for the mitochondrial chromosome used by the genome reference
# Commonly, either "chrM" or "MT"
mtChromosome="MT"

# ====== Peak calling ==========================================================
# Path to the peak calling directory
peakCallingDir="${projPath}/peak_calling"

# Effective genome size for the coverage files, used for normalization
# Refer to utils/resources
genomeSize="XXXXXXXXXX"   # <-- Update this

# q-value for MACS2 peak calling
qValue="0.01"

# p-value for MACS2 peak calling, not recommended to use
# Note that if p-value is set, q-value is ignored
pValue=""

# Whether to perform cutoff analysis, must be either "true" or "false"
# Note that, if true, the analysis tries to determine the cutoff value for the
# q-value, but it would require a significant higher amount of computational time
cutoffAnalysis="false"

# Whether to  build the shifting model during peak calling
# By default, a model is built by MACS2 to predict the shifting size
# If false, "--nomodel" is specified during peak calling, and "extSize" is used
useModel="true"

# Extension size for the MACS2 peak calling, the size extended to 3' of the read
# If "useModel" is false, this value is used
extsize=200

# ====== Processing called peaks ===============================================

#### Blacklist filtering

# Folder name of the raw peaks
# The peak calling step will generate a subfolder under the macs2 directory
# indicating the parameters used for peak calling
# Users must select and specify the folder name of the raw peaks that will be
# used for the downstream analysis
rawPeaks="q0.01_nolambda"  # <-- Update this

#### BAM to coverage

# Path to the coverage files
coverageFilePath="${projPath}/coverage_tracks"

# Coverage file format, either "bedgraph" or "bigwig"
outCoverageFormat="bigwig"

# Normalization method for the coverage files, one of RPKM, CPM, BPM, RPGC, None
coverageNorm="RPGC"

# Bin size for the coverage files
bamCoverageBinSize="20"

#### Compute matrix

# Reference point for the computeMatrix
referencePoint="center"

# Length of the region before the start of the peak
beforeRegionStartLength="3000"

# Length of the region after the end of the peak
afterRegionStartLength="3000"

# Length of the region before the start of the TSS
beforeTSSLength="1000"

# Length of the region after the end of the TSS
afterTSSLength="1000"




# [END] 