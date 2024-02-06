# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Trimmed fastq file directory: ${trimmedDir}"
echo "  bowtie2 Index: ${bowtieIndex}"
echo "  Alignment directory: ${alignmentDir}"
echo "  Maximum insert size: ${maxInsertLength}"
echo "  keep SAM files: ${keepSAM}"

# [Main]

# Align trimmed reads to the genome using Bowtie2
# We have assumed that the raw reads has been trimmed.
# Otherwise, replace "--local" with "--end-to-end"

# Check if the required directories exist, if not, create them
if [ ! -d ${alignmentDir} ]; then
    mkdir ${alignmentDir}
    mkdir ${alignmentDir}/sam
    mkdir ${alignmentDir}/sam_1000
    mkdir ${alignmentDir}/bam
    mkdir ${alignmentDir}/bowtie2_summary
fi

if [ ! -d ${alignmentDir}/sam ]; then
    mkdir ${alignmentDir}/sam
fi

if [ ! -d ${alignmentDir}/bam ]; then
    mkdir ${alignmentDir}/bam
fi

if [ ! -d ${alignmentDir}/bowtie2_summary ]; then
    mkdir ${alignmentDir}/bowtie2_summary
fi

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and align them using Bowtie2

for file in $(find ${trimmedDir} -type f -name '*_R1_*.fastq.gz'); do
    sample=$(basename $file _R1_trimmed.fastq.gz)

    echo "Aligning ${sample}..."
    echo "This may take a while..."

    # Step 1: Align trimmed reads to the genome using Bowtie2
    forward_file="${sample}_R1_trimmed.fastq.gz"
    reverse_file="${sample}_R2_trimmed.fastq.gz"
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 \
        -I 10 -X ${maxInsertLength} \
        -x ${bowtieIndex} \
        -p ${threads} \
        -1 ${trimmedDir}/${forward_file} \
        -2 ${trimmedDir}/${reverse_file} \
        -S ${alignmentDir}/sam/${sample}_bowtie2.sam \
          &> ${alignmentDir}/bowtie2_summary/${sample}_bowtie2.txt
    
    # Step 2: Sample top 1000 alignments
    head -n 1000 ${alignmentDir}/sam/${sample}_bowtie2.sam \
      > ${alignmentDir}/sam_1000/${sample}_bowtie2_top1000.sam
    
    # Step 3: Convert SAM to BAM
    samtools view -@ ${threads} -bS \
      ${alignmentDir}/sam/${sample}_bowtie2.sam \
        > ${alignmentDir}/bam/${sample}_bowtie2.bam
    
    # Step 4: Remove SAM file if not needed
    if [[ "${keepSAM}" =~ ^[f] ]]; then
        rm ${alignmentDir}/sam/${sample}_bowtie2.sam
    fi

    echo "Alignment of ${sample} complete"
done

# [END]