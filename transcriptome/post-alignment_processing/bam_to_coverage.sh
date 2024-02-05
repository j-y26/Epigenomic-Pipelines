# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "STAROutputDir: ${STAROutputDir}"
echo "coverageFilePath: ${coverageFilePath}"
echo "outCoverageFormat: ${outCoverageFormat}"
echo "coverageNorm: ${coverageNorm}"
echo "genomeSize: ${genomeSize}"
echo "binSize: ${binSize}"

# [Main]

# Create output directory if it doesn't exist

if [ ! -d ${coverageFilePath} ]; then
    mkdir ${coverageFilePath}
fi

# Output format suffix
if [ ${outCoverageFormat} = "bedgraph" ]; then
    outSuffix="bedgraph"
else
    outCoverageFormat="bigwig"
    outSuffix="bw"
fi

# Convert BAM to coverage file by iterating through each sample
# Do not extend reads for RNA-seq data

for file in $(find ${STAROutputDir} -type f -name '*.bam'); do
    sample=$(basename $file Aligned.sortedByCoord.out.bam)
    echo "Converting ${sample} BAM to coverage file"
    bamCoverage \
        -b ${STAROutputDir}/${sample}Aligned.sortedByCoord.out.bam \
        -o ${coverageFilePath}/${sample}.${outSuffix} \
        --outFileFormat ${outCoverageFormat} \
        --normalizeUsing ${coverageNorm} \
        --effectiveGenomeSize ${genomeSize} \
        --centerReads \
        --binSize ${binSize} \
        --numberOfProcessors ${threads}
    
# [END]