# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Filtered peak directory: ${peakCallingDir}/filtered_peaks/${rawPeaks}"
echo "  Coverage file path: ${coverageFilePath}"
echo "  GTF file: ${genomeGtfFile}"

# Create output directory if it doesn't exist
if [ ! -d ${peakCallingDir}/peak_coverage/${rawPeaks} ]; then
    mkdir -p ${peakCallingDir}/peak_coverage/${rawPeaks}
fi

if [ ! -d ${peakCallingDir}/tss_coverage/${rawPeaks} ]; then
    mkdir -p ${peakCallingDir}/tss_coverage/${rawPeaks}
fi

# Test whether any bigwig files exist in the coverage directory
bwFiles=$(find ${coverageFilePath} -name "*.bw")
if [ -z "${bwFiles}" ]; then
    echo "No bigwig files found in the coverage directory"
    echo "computeMatrix requires bigwig files associated with the peak files"
    exit 1
fi

# Obtain peak and coverage files
peakFiles=$(find ${peakCallingDir}/filtered_peaks/${rawPeaks} -name "*.bed" | sort)
labels=$(find ${peakCallingDir}/filtered_peaks/${rawPeaks} -name "*.bed" | sort | sed 's/.*\///' | sed 's/\.bed//')
bwFiles=$(find ${peakCallingDir}/filtered_peaks/${rawPeaks} -name "*.bed" | sort | sed 's|.*/\(.*\)\.bed|'"${coverageFilePath}"'/\1.bw|')

# Compute matrix for the peak files
echo "Computing matrix over peaks"
computeMatrix reference-point \
    --referencePoint ${referencePoint} \
    --regionsFileName ${peakFiles} \
    --scoreFileName ${bwFiles} \
    --outFileName ${peakCallingDir}/peak_coverage/${rawPeaks}/coverage_matrix_peaks.gz \
    --outFileNameMatrix ${peakCallingDir}/peak_coverage/${rawPeaks}/coverage_matrix_peaks.tab \
    --samplesLabel ${labels} \
    --beforeRegionStartLength ${beforeRegionStartLength} \
    --afterRegionStartLength ${afterRegionStartLength} \
    --skipZeros \
    --numberOfProcessors ${threads}
echo "Matrix computed"

# Compute matrix for TSS
echo "Computing matrix over TSS"
computeMatrix reference-point \
    --referencePoint TSS \
    --regionsFileName ${genomeGtfFile} \
    --scoreFileName ${bwFiles} \
    --outFileName ${peakCallingDir}/tss_coverage/${rawPeaks}/coverage_matrix_tss.gz \
    --outFileNameMatrix ${peakCallingDir}/tss_coverage/${rawPeaks}/coverage_matrix_tss.tab \
    --smartLabels \
    --beforeRegionStartLength ${beforeTSSLength} \
    --afterRegionStartLength ${afterTSSLength} \
    --skipZeros \
    --numberOfProcessors ${threads}
echo "Matrix computed"

# [END]