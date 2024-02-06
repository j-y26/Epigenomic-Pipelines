# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  alignmentDir: ${alignmentDir}"


# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${alignmentDir}/fragment_size ]; then
    mkdir ${alignmentDir}/fragment_size
fi


# Calculate the fragment size for each sample

for file in $(find ${alignmentDir}/bam_nodup -name "*.bam"); do
    sample=$(basename $file _nodups.bam)
    echo "Calculating fragment size for ${sample}"

    # Step 1: Extract the 9th column from the BAM file
    samtools view -F 0x04 -@ ${threads} ${alignmentDir}/bam_nodup/${sample}_nodups.bam | \

    # Step 2: Generate the absolute value of the fragment length
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \

    # Step 3: Sort the fragment length
    sort | \

    # Step 4: Count the number of occurrence of each fragment length
    uniq -c | \

    # Step 5: Print the fragment length and occurrence
    awk -v OFS="\t" '{print $2, $1/2}' > ${alignmentDir}/fragment_size/${sample}_fragment_size.txt
done

# [END]