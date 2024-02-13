# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  Mitochondrial reads QC directory: ${mt_reads}"


# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/mt_reads ]; then
    mkdir ${alignmentDir}/mt_reads
fi

mt_reads="${alignmentDir}/mt_reads"

# Check if the output csv file exists, if so, ask to delete it
if [ -f ${mt_reads}/mt_reads.csv ]; then
    echo "Existing ${mt_reads}/mt_reads.csv found, delete and proceed? (y/n)"
    read response
    if [ "$response" != "y" ]; then
        echo "Exiting..."
        exit 1
    fi
    echo "Proceeding..."
    rm ${mt_reads}/mt_reads.csv
    echo "Deleted existing ${mt_reads}/mt_reads.csv"
fi

# Initialize a CSV file for saving the number of mitochondrial reads
echo "sample,mt_reads,total_reads,mark" > ${mt_reads}/mt_reads.csv

# Iterate over all txt files (one for each mark) and perform analysis on all
# samples indicated in the file

for file in $(find ${markedSamples} -type f -name "samples_*.txt"); do
    mark=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')
    echo "Processing mark ${mark}"
    samples=$(cat ${file})
    bamFiles=""
    labels=""
    for sample in ${samples}; do
        bamFiles="${bamFiles} ${alignmentDir}/filtered_bam/${sample}${bamSuffix}"
        labels="${labels} ${sample}"
    done

    # Calculate the number of total and mitochondrial reads for each sample

    # Iterate over the index of bamFiles and labels
    for i in $(seq 0 $((${#bamFiles[@]}-1))); do
        sample=${labels[$i]}
        bamFile=${bamFiles[$i]}
        echo "Calculating mitochondrial/total reads for ${sample}"

        # Step 1: Extract the mitochondrial reads
        mt_reads=$(samtools view -c ${bamFile} 'chrM|chrMT|MT|chrMito|M')
        echo "  Mitochondrial reads: ${mt_reads} for ${sample}"

        # Step 2: Calculate the total reads
        total_reads=$(samtools view -c ${bamFile})
        echo "  Total reads: ${total_reads} for ${sample}"

        # Step 3: Append the results to the CSV file
        echo "${sample},${mt_reads},${total_reads},${mark}" >> ${mt_reads}/mt_reads.csv

        echo "Done for ${sample}"
    done
    echo "Done with mark ${mark}"
    echo ""
done

# [END]