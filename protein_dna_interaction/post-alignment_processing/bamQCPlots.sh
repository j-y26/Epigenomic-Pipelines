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
echo "  Correlation method: ${correlationMethod}"
echo "  Color map: ${colorMap}"
echo "  Correlation plot height: ${correlationPlotHeight}"
echo "  Correlation plot width: ${correlationPlotWidth}"
echo "  PCA plot height: ${pcaPlotHeight}"
echo "  PCA plot width: ${pcaPlotWidth}"
echo "  Number of top components for PCA: ${ntop}"

# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/bam_qc ]; then
    mkdir ${alignmentDir}/bam_qc
fi

# Coverage file calculated by multiBamSummary
coverageFile=${alignmentDir}/bam_qc/multiBamSummary_${binSize}.npz

# Correlation plot and PCA plot

# Iterate over all txt files (one for each mark) and perform analysis on all
# samples indicated in the file

for file in $(find ${markedSamples} -type f -name "samples_*.txt"); do
    mark=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')
    echo "Processing mark ${mark}"
    samples=$(cat ${file})
    labels=""
    for sample in ${samples}; do
        labels="${labels} ${sample}"
    done

    # Find the coverage file for the current mark
    coverageFile=${alignmentDir}/bam_qc/multiBamSummary_${mark}_${binSize}.npz

    # Check if the file exists, if not, skip the current mark
    if [ ! -f ${coverageFile} ]; then
        echo "Coverage file for mark ${mark} not found, skipping"
        continue
    fi

    # Correlation plot
    plotCorrelation \
        --corData ${coverageFile} \
        --corMethod ${correlationMethod} \
        --whatToPlot heatmap \
        --plotFile ${alignmentDir}/bam_qc/correlation_heatmap_${mark}.pdf \
        --outFileCorMatrix ${alignmentDir}/bam_qc/correlation_matrix_${mark}.tab \
        --skipZeros \
        --removeOutliers \
        --colorMap ${colorMap} \
        --plotNumbers \
        --labels ${labels} \
        --plotHeight ${correlationPlotHeight} \
        --plotWidth ${correlationPlotWidth}

    # PCA plot
    plotPCA \
        --corData ${coverageFile} \
        --plotFile ${alignmentDir}/bam_qc/pca_plot_${mark}.pdf \
        --labels ${labels} \
        --outFileNameData ${alignmentDir}/bam_qc/pca_data_${mark}.tab \
        --ntop ${ntop} \
        --plotHeight ${pcaPlotHeight} \
        --plotWidth ${pcaPlotWidth}
done

# [End]