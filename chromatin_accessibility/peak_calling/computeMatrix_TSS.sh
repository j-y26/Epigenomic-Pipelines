# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Coverage file path: ${coverageFilePath}"
echo "  GTF file: ${genomeGtfFile}"

# Create output directory if it doesn't exist
if [ ! -d ${peakCallingDir}/tss_coverage ]; then
    mkdir ${peakCallingDir}/tss_coverage
fi

# Test whether any bigwig files exist in the coverage directory
bwFiles=$(find ${coverageFilePath} -name "*.bw")
if [ -z "${bwFiles}" ]; then
    echo "No bigwig files found in the coverage directory"
    echo "computeMatrix requires bigwig files associated with the peak files"
    exit 1
fi

# Obtain peak and coverage files
bwFiles=$(find ${coverageFilePath} -name "*.bed" | sort | sed 's|.*/\(.*\)\.bed|'"${coverageFilePath}"'/\1.bw|')

# Compute matrix for the peak files
computeMatrix reference-point \
    --referencePoint TSS \
    --regionsFileName ${genomeGtfFile} \
    --scoreFileName ${bwFiles} \
    --outFileName ${peakCallingDir}/tss_coverage/coverage_matrix_tss.gz \
    --outFileNameMatrix ${peakCallingDir}/tss_coverage/coverage_matrix_tss.tab \
    --smartLabels \
    --beforeRegionStartLength ${beforeTSSLength} \
    --afterRegionStartLength ${afterTSSLength} \
    --skipZeros \
    --numberOfProcessors ${threads}

# [END]