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

# Check whether rRNA cleaning is enabled
if [ ${rRNAFiltering} == "true" ]; then
    trimmedDir=${rRNAcleanDir}
    suffix="cleaned.fastq.gz"
else
    suffix="trimmed.fastq.gz"
fi

# Align the reads to the reference genome
# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# We use trimmed files for alignment, so we need to loop through the trimmed directory

for file in $(find ${trimmedDir} -type f -name '*R1_*.fastq.gz'); do
    sample=$(basename $file _R1_${suffix})
    echo "Aligning ${sample} to reference genome..."
    echo "This may take a while..."
    forward_file="${sample}_R1_${suffix}"
    reverse_file="${sample}_R2_${suffix}"

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
    --readFilesIn ${trimmedDir}/${forward_file} ${trimmedDir}/${reverse_file}
    echo "Alignment for ${sample} complete"
done

# [END]
