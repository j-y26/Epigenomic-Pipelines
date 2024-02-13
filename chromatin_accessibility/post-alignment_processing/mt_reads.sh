# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  Mitochondrial reads QC directory: ${alignmentDir}/mt_reads"


# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/mt_reads ]; then
    mkdir ${alignmentDir}/mt_reads
fi

# Check if the output csv file exists, if so, ask to delete it
if [ -f ${alignmentDir}/mt_reads/mt_reads.csv ]; then
    echo "Existing ${alignmentDir}/mt_reads/mt_reads.csv found, delete and proceed? (y/n)"
    read response
    if [ "$response" != "y" ]; then
        echo "Exiting..."
        exit 1
    fi
    echo "Proceeding..."
    rm ${alignmentDir}/mt_reads/mt_reads.csv
    echo "Deleted existing ${alignmentDir}/mt_reads/mt_reads.csv"
fi

# Initialize a CSV file for saving the number of mitochondrial reads
echo "sample,mt_reads,total_reads" > ${alignmentDir}/mt_reads/mt_reads.csv

# Calculate the number of total and mitochondrial reads for each sample

for file in $(find ${alignmentDir}/filtered_bam -name "*.bam"); do
    sample=$(basename $file .bam)
    echo "Calculating mitochondrial/total reads for ${sample}"

    # Step 1: Extract the mitochondrial reads
    mt_reads=$(samtools view -c ${file} ${mtChromosome})
    echo "  Mitochondrial reads: ${mt_reads} for ${sample}"

    # Step 2: Calculate the total reads
    total_reads=$(samtools view -c ${file})
    echo "  Total reads: ${total_reads} for ${sample}"

    # Step 3: Append the results to the CSV file
    echo "${sample},${mt_reads},${total_reads}" >> ${alignmentDir}/mt_reads/mt_reads.csv

    echo "Done for ${sample}"
done

# [END]