# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  PICARD_JAR: ${PICARD_JAR}"

# [Main]

# Check if the alignment directory exists, if not, create it
if [ ! -d ${alignmentDir}/bam_nodup ]; then
    mkdir ${alignmentDir}/bam_nodup
fi


# Remove duplicates from bowtie2 aligned bam files

for file in $(find ${alignmentDir}/bam -type f -name '*_bowtie2.bam'); do
    sample=$(basename $file _bowtie2.bam)
    echo "Processing sample ${sample}..."

    # Step 1: Sort the bam file by coordinates
    echo "Sorting BAM file for ${sample}..."
    samtools sort -@ $threads -o ${alignmentDir}/bam/${sample}_coord_sorted.bam \
        ${alignmentDir}/bam/${sample}_bowtie2.bam
    echo "BAM file sorted"

    # Step 2: Mark duplicates
    echo "Marking duplicates for ${sample}..."
    java -jar ${PICARD_JAR} MarkDuplicates \
        I=${alignmentDir}/bam/${sample}_coord_sorted.bam \
        O=${alignmentDir}/bam_nodup/${sample}_markdups.bam \
        M=${alignmentDir}/bam_nodup/${sample}_dups_metrics.txt \
        REMOVE_DUPLICATES=false
    echo "Duplicates marked"

    # Step 3: Remove duplicates
    samtools view -@ $threads -b -F 0x0400 \
        ${alignmentDir}/bam_nodup/${sample}_markdups.bam > \
        ${alignmentDir}/bam_nodup/${sample}_nodup.bam

    # Step 4: Remove intermediate files
    echo "Removing intermediate files for ${sample}..."
    rm ${alignmentDir}/bam/${sample}_coord_sorted.bam
    rm ${alignmentDir}/bam/${sample}_markdups.bam
    echo "Intermediate files removed"

    echo "Sample ${sample} deduplication complete"
done

# [END]