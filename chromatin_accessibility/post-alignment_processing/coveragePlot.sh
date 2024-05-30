# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  Number of 1bp regions to sample: ${numberOfSamples}"
echo "  Blacklist file: ${blacklistFile}"

# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/bam_qc ]; then
    mkdir ${alignmentDir}/bam_qc
fi

# BAM files and labels
bamFiles=$(find ${alignmentDir}/filtered_bam -name "*.bam" | sort)
labels=$(find ${alignmentDir}/filtered_bam -name "*.bam" | sort | sed 's/.*\///' | sed 's/\.bam//')

# Generate coverage plot

plotCoverage \
    --bamfiles ${bamFiles} \
    --labels ${labels} \
    --plotFile ${alignmentDir}/bam_qc/coverage_plot.pdf \
    --outRawCounts ${alignmentDir}/bam_qc/coverage_counts.tab \
    --outCoverageMetrics ${alignmentDir}/bam_qc/coverage_metrics.tab \
    --skipZeros \
    --minMappingQuality 20 \
    --plotFileFormat pdf \
    --numberOfSamples ${numberOfSamples} \
    --blackListFileName ${blacklistFile} \
    --plotHeight ${coveragePlotHeight} \
    --plotWidth ${coveragePlotWidth} \
    --numberOfProcessors ${threads}

# [END]