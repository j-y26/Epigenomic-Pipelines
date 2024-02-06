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
projPath="/your/project/path"

# Path to the directory containing the FASTQ files
fastqDir="${projPath}/fastq" 

# Lane identifier
lane="L001" 

# Number of threads to use for the analysis (check your CPU for this value)
threads="8" 

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
genomeFastaFile="/path/to/genome.fa"

# Path to the genome GTF file
genomeGtfFile="/path/to/genome.gtf"

# Path to the output directory for the alignment
alignmentDir="${projPath}/alignment"

# Maximum insert length for paired-end reads
maxInsertLength="1000"

# Whether to preserve the human-readable SAM files, must be their "true" or "false"
keepSAM="false"

# ====== Post-alignment processing =============================================
# Path to Picard JAR file
PICARD_JAR="/path/to/picard.jar"





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