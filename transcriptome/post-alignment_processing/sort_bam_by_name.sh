# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "STAROutputDir: ${STAROutputDir}"
echo "nameSortedDir: ${nameSortedDir}"

# [Main]

# Create output directory if it doesn't exist
if [ ! -d ${nameSortedDir} ]; then
    mkdir ${nameSortedDir}
fi

# Sort the BAM files by name using samtools

for file in $(find ${STAROutputDir} -type f -name '*.bam'); do
    sample=$(basename $file Aligned.sortedByCoord.out.bam)
    echo "Sorting ${sample} by name..."
    samtools sort -n -@ ${threads} -O bam \
    -o ${nameSortedDir}/${sample}_name_sorted.bam ${file}
    echo "Sorting for ${sample} complete"
done

# [END] 