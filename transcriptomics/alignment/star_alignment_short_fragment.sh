# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "STARIndexDir: ${STARIndexDir}"
echo "STAROutputDir: ${STAROutputDir}"

echo "Configuring STAR for short fragments..."
echo "Fragment Length Filtering Parameters:"
echo "outFilterScoreMinOverLread: ${filterLreadScore}"
echo "outFilterMatchNminOverLread: ${filterLreadMatch}"
echo "outFilterMismatchNmax: ${filterMismatchNmax}"

if [ ${useTempFIFO} == "true" ]; then
    echo "Using temporary FIFO directory: ${tempFIFODir}"

else
    echo "Using STAROutputDir as temporary FIFO directory"
fi

# [Main]

# Check to ensure that the reference genome directory exists
if [ ! -d ${STARIndexDir} ]; then
    echo "Error: STAR index directory does not exist." >&2
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d ${STAROutputDir} ]; then
    mkdir ${STAROutputDir}
fi

# Align the reads to the reference genome
# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# We use trimmed files for alignment, so we need to loop through the trimmed directory

for file in $(find ${trimmedDir} -type f -name '*R1_trimmed.fastq.gz'); do
    sample=$(basename $file _R1_trimmed.fastq.gz)
    echo "Aligning ${sample} to reference genome..."
    echo "This may take a while..."
    forward_file="${sample}_R1_trimmed.fastq.gz"
    reverse_file="${sample}_R2_trimmed.fastq.gz"

    # Create a temporary directory for the FIFO files
    if [ ${useTempFIFO} == "true" ]; then
        outTempDir="--outTmpDir ${tempFIFODir}/${sample}"
    else
        outTempDir=""
    fi

    STAR \
    --runThreadN ${threads} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${STAROutputDir}/${sample} \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir ${STARIndexDir} \
    ${outTempDir} \
    --readFilesIn ${trimmedDir}/${forward_file} ${trimmedDir}/${reverse_file} \
    --outFilterScoreMinOverLread ${filterLreadScore} \
    --outFilterMatchNminOverLread ${filterLreadMatch} \
    --outFilterMismatchNmax ${filterMismatchNmax} \
    echo "Alignment for ${sample} complete"

    # Remove the temporary directory if it was created
    if [ ${useTempFIFO} == "true" ]; then
        rm -r ${tempFIFODir}/${sample}
    fi
done

# [END]
