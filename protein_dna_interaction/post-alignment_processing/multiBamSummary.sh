# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  binSize: ${binSize}"

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

    # Run multiBamSummary by deepTools on all samples
    # Operate on bins for the entire genome
    multiBamSummary bins \
        --bamfiles ${bamFiles} \
        --labels ${labels} \
        --numberOfProcessors ${threads} \
        --binSize ${binSize} \
        --outFileName ${alignmentDir}/bam_qc/multiBamSummary_${mark}_${binSize}.npz
done

# [END]
