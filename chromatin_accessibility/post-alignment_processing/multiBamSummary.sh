# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  Bin size: ${binSize}"

# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/bam_qc ]; then
    mkdir ${alignmentDir}/bam_qc
fi

# BAM files and labels
bamFiles=$(find ${alignmentDir}/filtered_bam -name "*.bam" | sort)
labels=$(find ${alignmentDir}/filtered_bam -name "*.bam" | sort | sed 's/.*\///' | sed 's/\.bam//')

# Run multiBamSummary by deepTools on all samples
# Operate on bins for the entire genome

multiBamSummary bins \
    --bamfiles ${bamFiles} \
    --labels ${labels} \
    --numberOfProcessors ${threads} \
    --binSize ${binSize} \
    --verbose \
    --outFileName ${alignmentDir}/bam_qc/multiBamSummary_${binSize}.npz

# [END]
