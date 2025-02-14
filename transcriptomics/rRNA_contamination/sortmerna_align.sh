# # [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "rRNAFastaFile: ${rRNAFastaFile}"
echo "sortmernaRef: ${sortmernaRef}"
echo "rRNAalignDir: ${rRNAalignDir}"
echo "rRNAcleanDir: ${rRNAcleanDir}"
echo "sortmernaWorkDir: ${sortmernaWorkDir}"

# [Main]

# Create output directory if it doesn't exist
if [ ! -d ${rRNAalignDir} ]; then
    mkdir -p ${rRNAalignDir}
fi

if [ ! -d ${rRNAcleanDir} ]; then
    mkdir -p ${rRNAcleanDir}
fi

# Check if the sortmernaWorkDir exists
if [ ! -d ${sortmernaWorkDir} ]; then
    mkdir -p ${sortmernaWorkDir}
fi

# Check if the reference files exist
if [ ! -f ${rRNAFastaFile} ]; then
    echo "Genome specific rRNA FASTA file not found: ${rRNAFastaFile}"
    exit 1
fi

if [ ! -f ${sortmernaRef} ]; then
    echo "SortMeRNA reference file not found: ${sortmernaRef}"
    exit 1
fi

# Check if key-value database already exists, if so, purge it
if [ -d ${sortmernaWorkDir}/kvdb ]; then
    rm -rf ${sortmernaWorkDir}/kvdb
fi

# Check if pigz is installed
if ! command -v pigz &> /dev/null; then
    echo "pigz could not be found"
    echo "Please install pigz before running this script"
    exit 1
fi

# Align the reads to the rRNA database
for file in $(find ${trimmedDir} -type f -name '*R1_trimmed.fastq.gz'); do
    sample=$(basename $file _R1_trimmed.fastq.gz)
    echo "Aligning ${sample} to rRNA reference..."
    echo "This may take a while..."
    forward_file="${sample}_R1_trimmed.fastq.gz"
    reverse_file="${sample}_R2_trimmed.fastq.gz"

    # Clean up the kvdb and readb directiories
    if [ -d ${sortmernaWorkDir}/kvdb ]; then
        rm -rf ${sortmernaWorkDir}/kvdb
    fi
    if [ -d ${sortmernaWorkDir}/readb ]; then
        rm -rf ${sortmernaWorkDir}/readb
    fi

    # Run SortMeRNA
    sortmerna \
      --ref ${rRNAFastaFile} \
      --ref ${sortmernaRef} \
      --reads ${trimmedDir}/${forward_file} \
      --reads ${trimmedDir}/${reverse_file} \
      --aligned ${rRNAalignDir}/${sample}_aligned \
      --other ${rRNAcleanDir}/${sample}_cleaned \
      --fastx \
      --out2 \
      --sam \
      --SQ \
      --num_alignments 1 \
      --workdir ${sortmernaWorkDir} \
      --threads ${threads}

    # Clean up the kvdb and readb directiories
    rm -rf ${sortmernaWorkDir}/kvdb
    rm -rf ${sortmernaWorkDir}/readb

    echo "SortMeRNA alignment for ${sample} complete"

    # Rename the output files to follow standard naming convention
    mv ${rRNAalignDir}/${sample}_aligned_fwd.fq.gz ${rRNAalignDir}/${sample}_R1_aligned.fastq.gz
    mv ${rRNAalignDir}/${sample}_aligned_rev.fq.gz ${rRNAalignDir}/${sample}_R2_aligned.fastq.gz

    mv ${rRNAcleanDir}/${sample}_cleaned_fwd.fq.gz ${rRNAcleanDir}/${sample}_R1_cleaned.fastq.gz
    mv ${rRNAcleanDir}/${sample}_cleaned_rev.fq.gz ${rRNAcleanDir}/${sample}_R2_cleaned.fastq.gz
done

# [END]