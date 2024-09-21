# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "genomeFastaFile: ${genomeFastaFile}"
echo "genomeGtfFile: ${genomeGtfFile}"
echo "STARIndexDir: ${STARIndexDir}"

# [Main]

# Build the index required for STAR alignment

# Create output directory if it doesn't exist
if [ ! -d ${STARIndexDir} ]; then
    mkdir -p ${STARIndexDir}
fi

# Create STAR index for the aligner

STAR --runMode genomeGenerate \
        --runThreadN ${threads} \
        --genomeDir ${STARIndexDir} \
        --genomeFastaFiles ${genomeFastaFile} \
        --sjdbGTFfile ${genomeGtfFile} \
        --sjdbOverhang 100

# [END]