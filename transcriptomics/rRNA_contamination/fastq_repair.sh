# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  rRNAcleanDir: ${rRNAcleanDir}"
echo "  Repair reports will be saved to ${rRNAcleanDir}/reports"

# [Main]

# Repairing rRNA-filtered fastq files with BBMap's repair.sh

# Create output directory if it doesn't exist

if [ ! -d ${rRNAcleanDir}/singletons ]; then
    mkdir -p ${rRNAcleanDir}/singletons
fi

if [ ! -d ${rRNAcleanDir}/reports ]; then
    mkdir -p ${rRNAcleanDir}/reports
fi

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and perform repairing

for file in $(find ${rRNAcleanDir} -type f -name '*_R1_cleaned.fastq.gz'); do
    sample=$(basename $file _R1_cleaned.fastq.gz)
    echo "Repairing ${sample}"
    forward_file="${sample}_R1_cleaned.fastq.gz"
    reverse_file="${sample}_R2_cleaned.fastq.gz"
    
    # Repair with BBMap's repair.sh
    repair.sh \
        in1=${rRNAcleanDir}/${forward_file} \
        in2=${rRNAcleanDir}/${reverse_file} \
        out1=${rRNAcleanDir}/${sample}_R1_cleaned_repaired.fastq.gz \
        out2=${rRNAcleanDir}/${sample}_R2_cleaned_repaired.fastq.gz \
        outs=${rRNAcleanDir}/singletons/${sample}_cleaned_singletons.fastq.gz 2>&1 | \
        tee ${rRNAcleanDir}/reports/${sample}_repair_report.txt

    echo "Done repairing ${sample}"
done

# [END]