# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Filtered peak directory: ${peakCallingDir}/filtered_peaks"
echo "  Coverage file path: ${coverageFilePath}"
echo "  GTF file: ${genomeGtfFile}"

# Create output directory if it doesn't exist
if [ ! -d ${peakCallingDir}/peak_coverage ]; then
    mkdir ${peakCallingDir}/peak_coverage
fi

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

# Compute matrix over peaks for each mark
# Here, we assume that the peak and coverage files are named by
# label.bed and label.bw, respectively
for file in $(find ${markedSamples} -type f -name "samples_*.txt"); do
    mark=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')
    echo "Processing mark ${mark}"
    samples=$(cat ${file})
    peakFiles=""
    labels=""
    bwFiles=""
    for sample in ${samples}; do
        peakFiles="${peakFiles} $(find ${peakCallingDir}/filtered_peaks -name "${sample}.bed")"
        labels="${labels} ${sample}"
        bwFiles="${bwFiles} $(find ${coverageFilePath} -name "${sample}.bw")"
    done
    peakFiles=$(echo ${peakFiles} | tr ' ' '\n' | sort | tr '\n' ' ')
    labels=$(echo ${labels} | tr ' ' '\n' | sort | tr '\n' ' ')
    bwFiles=$(echo ${bwFiles} | tr ' ' '\n' | sort | tr '\n' ' ')

    # Compute matrix for the peak files
    echo "Computing matrix over peaks"
    computeMatrix reference-point \
        --referencePoint ${referencePoint} \
        --regionsFileName ${peakFiles} \
        --scoreFileName ${bwFiles} \
        --outFileName ${peakCallingDir}/peak_coverage/${mark}_coverage_matrix_peaks.gz \
        --outFileNameMatrix ${peakCallingDir}/peak_coverage/${mark}_coverage_matrix_peaks.tab \
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
        --outFileName ${peakCallingDir}/tss_coverage/${mark}_coverage_matrix_tss.gz \
        --outFileNameMatrix ${peakCallingDir}/tss_coverage/${mark}_coverage_matrix_tss.tab \
        --smartLabels \
        --beforeRegionStartLength ${beforeTSSLength} \
        --afterRegionStartLength ${afterTSSLength} \
        --skipZeros \
        --numberOfProcessors ${threads}
    echo "Matrix computed"

done

# [END]