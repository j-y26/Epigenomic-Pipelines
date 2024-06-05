# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Input peak directory: ${peakCallingDir}/macs2/${rawPeaks}"
echo "  Filtered peak output directory: ${peakCallingDir}/filtered_peaks"
echo "  Blacklist file: ${blacklistFile}"


# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${peakCallingDir}/filtered_peaks ]; then
    mkdir ${peakCallingDir}/filtered_peaks
fi

# Filter the peaks using the blacklist file
echo ""
echo "Filtering peaks using the blacklist file"
echo "----------------------------------------"

# Narrow peak files
for file in $(find ${peakCallingDir}/macs2/${rawPeaks} -type f -name '*.narrowPeak'); do
    sample=$(basename $file _peaks.narrowPeak)
    echo "Filtering peaks for ${sample}..."

    # Filter the peaks
    bedtools intersect -v -a ${file} -b ${blacklistFile} > \
      ${peakCallingDir}/filtered_peaks/${sample}.bed
    
    echo "Peaks filtered for ${sample}"
done

# Broad peak files
for file in $(find ${peakCallingDir}/macs2/${rawPeaks} -type f -name '*.broadPeak'); do
    sample=$(basename $file _peaks.broadPeak)
    echo "Filtering peaks for ${sample}..."

    # Filter the peaks
    bedtools intersect -v -a ${file} -b ${blacklistFile} > \
      ${peakCallingDir}/filtered_peaks/${sample}.bed
    
    echo "Peaks filtered for ${sample}"
done

echo "Peak filtering complete"

# [END]