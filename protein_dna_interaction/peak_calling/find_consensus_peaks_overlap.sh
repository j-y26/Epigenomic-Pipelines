# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Filtered peak directory: ${peakCallingDir}/filtered_peaks/${rawPeaks}"
echo "  Consensus peak directory: ${peakCallingDir}/consensus_peaks/${rawPeaks}"
echo "  Minimum fraction of peaks needs to be overlapped: ${minOverlapFraction}"

# Create output directory if it doesn't exist
if [ ! -d ${peakCallingDir}/consensus_peaks/${rawPeaks} ]; then
    mkdir -p ${peakCallingDir}/consensus_peaks/${rawPeaks}
fi

# Iterate over the groups of samples identifying the replicate names
for file in $(find ${peakCallingDir}/consensus_peaks -type f -name "samples_*.txt"); do
    group=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')

    echo "Finding consensus peaks for ${group}"
    samples=$(cat ${file})
    peakFiles=""
    for sample in ${samples}; do
        peakFiles="${peakFiles} $(find ${peakCallingDir}/filtered_peaks/${rawPeaks} -name "${sample}.bed")"
    done
    peakFiles=$(echo ${peakFiles} | tr ' ' '\n' | sort | tr '\n' ' ')

    # Find consensus peaksgroup=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')

    # Get the first peak file name in the list
    firstPeakFile=$(echo ${peakFiles} | cut -d ' ' -f 1)
    # The rest of the peak files
    restPeakFiles=$(echo ${peakFiles} | cut -d ' ' -f 2-)

    # Find the consensus peaks
    bedtools intersect -a ${firstPeakFile} -b ${restPeakFiles} -f ${minOverlapFraction} -r -wa -u > \
        ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}_consensus_overlap.bed

done

echo "Consensus peaks have been found and saved in ${peakCallingDir}/consensus_peaks/${rawPeaks}"

# [END]