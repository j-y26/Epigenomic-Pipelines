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

# Global parameters
projPath="/your/project/path" # Path to the project home directory
fastqDir="${projPath}/fastq" # Path to the directory containing the FASTQ files
lane="L001" # Lane identifier
threads="8" # Number of threads to use for the analysis

# Data preprocessing
fastqcOutDir="${projPath}/fastqc" # Path to the directory containing the FASTQC output
TRIMMOMATIC_JAR="/path/to/trimmomatic.jar" # Path to the Trimmomatic JAR file
TRIMMOMATIC_ADAPTER_FILE="/path/to/adapter.fa" # Path to the Trimmomatic adapter file
# Note that the adapter must be chosen corresponding to the library generation method
# These fasta files are usually stored along with the Trimmomatic JAR file

trimmedDir="${projPath}/fastq_trimmed" # Path to the directory containing the trimmed FASTQ files

