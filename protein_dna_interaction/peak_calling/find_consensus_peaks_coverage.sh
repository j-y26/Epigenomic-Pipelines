# Finding consensus peaks from multiple peak calling results
# Consensus peaks are found by calling peaks over overlapping regions of minimum
# number of replicates.

# The script is adapted from taoliu/merge_then_call_consensus.sh

# This script will find the consensus peak regions from peak files (in
# BED format) of multiple samples by:
#
# 1. Converting the peak file of each sample into non-overlapping 3
# cols BED file and concatenating them;
#
# 2. Sorting the concatenated file and Building a genome coverage
# track in BedGraph, of which the value (the 3rd col) indicates the
# number of samples with a peak at a particular location;
#
# 3. Using MACS to call regions being covered by peaks from more than
# a certain number of samples.


# ------------------------------------------------------------------

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
echo "  Genome size file: ${genomeSizeFile}"
echo "  Minimum number of replicates identified as consensus peaks: ${minReplicates}"
echo "  Minimum peak size: ${minPeakSize}"
echo "  Maximum gap size: ${maxGapSize}"

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

    # Find consensus peaks

    # 1. Concatenate the non-overlapping peak files of replicates
    for f in ${peakFiles}; do
        bedtools sort -i ${f} | bedtools merge -i - | cut -f 1,2,3 >> \
          ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}.all.bed
    done

    # 2. Sort the concatenated file and build a genome coverage track, of which
    # the value (the 3rd col) indicates the number of replicates with a peak at
    # a particular location
    bedtools sort -i ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}.all.bed | \
        bedtools genomecov -bga -i - -g ${genomeSizeFile} > \
        ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}_coverage.bdg
    
    # 3. Call consensus peaks
    macs2 bdgpeakcall -i ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}_coverage.bdg \
        -o ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}_consensus_cov.bed \
        --no-trackline \
        --cutoff ${minReplicates} --min-len ${minPeakSize} --max-gap ${maxGapSize}

    # Remove the intermediate files
    rm -f ${peakCallingDir}/consensus_peaks/${rawPeaks}/${group}.all.bed

done

echo "Consensus peaks have been found and saved in ${peakCallingDir}/consensus_peaks/${rawPeaks}"

# [END]