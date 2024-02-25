# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Filtered BAM file directory: ${alignmentDir}/filtered_bam"
echo "  Coverage file path: ${coverageFilePath}"
echo "  Coverage file format: ${outCoverageFormat}"
echo "  Normalization method: ${coverageNorm}"
echo "  Bin size: ${bamCoverageBinSize}"
echo "  Genome size: ${genomeSize}"
echo "  Threads: ${threads}"


# [Main]

# Create output directory if it doesn't exist

if [ ! -d ${coverageFilePath} ]; then
    mkdir ${coverageFilePath}
fi

# Output format suffix
if [ ${outCoverageFormat} = "bedgraph" ]; then
    outSuffix="bdg"
else
    outCoverageFormat="bigwig"
    outSuffix="bw"
fi

# Convert BAM to coverage file by iterating through each sample

for file in $(find ${alignmentDir}/filtered_bam -type f -name '*.bam'); do
    sample=$(basename $file .bam)
    echo "Converting ${sample} BAM to coverage file"
    bamCoverage \
        -b ${file} \
        -o ${coverageFilePath}/${sample}.${outSuffix} \
        --outFileFormat ${outCoverageFormat} \
        --normalizeUsing ${coverageNorm} \
        --effectiveGenomeSize ${genomeSize} \
        --extendReads \
        --ignoreDuplicates \
        --centerReads \
        --binSize ${bamCoverageBinSize} \
        --numberOfProcessors ${threads}
    
# [END]