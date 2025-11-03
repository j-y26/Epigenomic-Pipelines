# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"

# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/filtered_bam ]; then
    mkdir ${alignmentDir}/filtered_bam
fi


# Execute the following steps in order for each bam file:
# 1. Filter out unmapped reads and keep only properly paired reads
# 2. Sort the reads by coordinate
# 3. Index the sorted reads

for file in $(find ${alignmentDir}/bam_nodup -type f -name "*.bam"); do
    sample=$(basename $file _nodup.bam)
    echo "Filtering and indexing ${sample}..."

    # Step 1: Filter out unmapped reads and keep only properly paired reads
    echo "  Filtering out unmapped reads and keeping only properly paired reads..."
    samtools view -@ ${threads} -F 0x04 -f 0x02 -b ${file} > ${alignmentDir}/filtered_bam/${sample}_filtered.bam
    echo "  Done."

    # Step 2: Sort the reads by coordinate
    echo "  Sorting the reads by coordinate..."
    samtools sort -@ ${threads} ${alignmentDir}/filtered_bam/${sample}_filtered.bam -o ${alignmentDir}/filtered_bam/${sample}.bam
    echo "  Done."

    # Step 3: Index the sorted reads
    echo "  Indexing the sorted reads..."
    samtools index -@ ${threads} ${alignmentDir}/filtered_bam/${sample}.bam
    echo "  Done."

    # Remove intermediate files
    rm ${alignmentDir}/filtered_bam/${sample}_filtered.bam

    echo "Processing ${sample} complete." 
done

# [END]