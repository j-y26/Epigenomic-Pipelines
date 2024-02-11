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
# Iterate over all txt files (one for each mark) and perform analysis on all
# samples indicated in the file

for file in $(find ${markedSamples} -type f -name "samples_*.txt"); do
    mark=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')
    echo "Processing mark ${mark}"
    samples=$(cat ${file})
    bamFiles=""
    labels=""
    for sample in ${samples}; do
        bamFiles="${bamFiles} ${alignmentDir}/bam_nodup/${sample}${bamSuffix}"
        labels="${labels} ${sample}"
    done

    # Calculate the fragment size distribution
    # Specifically, reads that overlap with the blacklist regions are removed
    bamPEFragmentSize \
        --bamfiles ${bamFiles} \
        --samplesLabel ${labels} \
        --histogram ${alignmentDir}/fragment_size/frag_dist_${mark}_deeptools.pdf \
        --plotFileFormat pdf \
        --outRawFragmentLengths ${alignmentDir}/fragment_size/frag_len_${mark}_deeptools.txt \
        --maxFragmentLength ${maxInsertLength} \
        --numberOfProcessors ${threads} \
        --blackListFileName ${blacklistFile} \
        --verbose
done
echo "Fragment size analysis complete"

# [End]