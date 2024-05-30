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

# Iterate over all txt files (one for each mark) and perform analysis on all
# samples indicated in the file

for file in $(find ${markedSamples} -type f -name "samples_*.txt"); do
    mark=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')
    echo "Processing mark ${mark}"
    samples=$(cat ${file})
    bamFiles=""
    labels=""
    for sample in ${samples}; do
        bamFiles="${bamFiles} ${alignmentDir}/filtered_bam/${sample}.bam"
        labels="${labels} ${sample}"
    done

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
        --numberOfProcessors ${threads} \
        --verbose
done

# [END]