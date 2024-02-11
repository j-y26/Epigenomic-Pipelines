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

# Sample labels
labels=$(find ${alignmentDir}/filtered_bam -name "*.bam" | sort | sed 's/.*\///' | sed 's/\.bam//')

# Correlation plot and PCA plot

# Correlation plot
plotCorrelation \
    --corData ${coverageFile} \
    --corMethod ${correlationMethod} \
    --whatToPlot heatmap \
    --plotFile ${alignmentDir}/bam_qc/correlation_heatmap.pdf \
    --outFileCorMatrix ${alignmentDir}/bam_qc/correlation_matrix.tab \
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
    --plotFile ${alignmentDir}/bam_qc/pca_plot.pdf \
    --labels ${labels} \
    --outFileNameData ${alignmentDir}/bam_qc/pca_data.tab \
    --ntop ${ntop} \
    --plotHeight ${pcaPlotHeight} \
    --plotWidth ${pcaPlotWidth}

# [End]