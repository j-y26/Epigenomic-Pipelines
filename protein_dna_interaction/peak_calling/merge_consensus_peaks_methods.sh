# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Consensus peak directory: ${peakCallingDir}/consensus_peaks/${rawPeaks}"

# Check if consensus peaks identified by both methods exist
if [ ! -d ${peakCallingDir}/consensus_peaks/${rawPeaks} ]; then
    echo "Consensus peaks have not been found for ${rawPeaks}"
    exit 1
fi

# Peaks from both methods identified by _consensus_cov.bed and _consensus_overlap.bed
covPeaks=$(find ${peakCallingDir}/consensus_peaks/${rawPeaks} -name "*_consensus_cov.bed")
overlapPeaks=$(find ${peakCallingDir}/consensus_peaks/${rawPeaks} -name "*_consensus_overlap.bed")

# Check if consensus peaks identified by both methods exist
if [ -z "${covPeaks}" ] || [ -z "${overlapPeaks}" ]; then
    echo "Consensus peaks have not been found for ${rawPeaks}"
    exit 1
fi

# Merge the consensus peaks identified by the two methods
for covPeak in ${covPeaks}; do
    group=$(basename ${covPeak} | sed 's/_consensus_cov.bed//g')
    overlapPeak=$(find ${peakCallingDir}/consensus_peaks/${rawPeaks} -name "${group}_consensus_overlap.bed")

    if [ -z "${overlapPeak}" ]; then
        echo "Consensus peaks have not been found for ${group}"
    else
        echo "Merging consensus peaks for ${group}"
        cat ${covPeak} ${overlapPeak} | \
            cut -f 1-3 | \
            bedtools sort -i - | \
            bedtools merge -i - > \
            ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}_consensus_merged.bed
    fi
done

# [END]