# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  Fragment size output directory: ${alignmentDir}/fragment_size"
echo "  Blacklist file: ${blacklistFile}"
echo "  Max insert length: ${maxInsertLength}"

# [Main]

# Check if the output directory exists, if not, create it
if [ ! -d ${alignmentDir}/fragment_size ]; then
    mkdir -p ${alignmentDir}/fragment_size
fi

# Using deepTools bamPEFragmentSize to estimate the fragment size of the paired-end reads

# Obtain all BAM files and respective labels, where labels are the sample names
# with the suffix .bam and the path to the BAM file removed
bamFiles=$(ls ${alignmentDir}/filtered_bam/*.bam | sort)
labels=$(ls ${alignmentDir}/filtered_bam/*.bam | sort | sed 's/.*\///g' | sed 's/.bam//g')

# Calculate the fragment size distribution
# Specifically, reads that overlap with the blacklist regions are removed

bamPEFragmentSize \
    --bamfiles ${bamFiles} \
    --samplesLabel ${labels} \
    --histogram ${alignmentDir}/fragment_size/frag_dist_deeptools.pdf \
    --plotFileFormat pdf \
    --outRawFragmentLengths ${alignmentDir}/fragment_size/frag_len_deeptools.txt \
    --maxFragmentLength ${maxInsertLength} \
    --numberOfProcessors ${threads} \
    --blackListFileName ${blacklistFile} \
    --verbose

# [End]