# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Alignment directory: ${alignmentDir}"
echo "  Peak calling directory: ${peakCallingDir}"

# [Main]

# Check if the directory exists, if not, create it
if [ ! -d ${peakCallingDir} ]; then
    mkdir ${peakCallingDir}
    mkdir ${peakCallingDir}/macs2
fi

if [ ! -d ${peakCallingDir}/macs2 ]; then
    mkdir ${peakCallingDir}/macs2
fi

# Whether to perform cutoff analysis
if [${cutoffAnalysis} == "true"]; then
    echo "Cutoff analysis is enabled"
    cutoff="--cutoff-analysis"
else
    cutoff=""
fi

# ATAC-seq peaks are cumulated as narrow regions in accessible promoters and 
# enhancers, so the narrow peak mode is used
# Since ATACseq does not have control samples, background lambda is used

# Iterate over all samples and perform peak calling with MACS2

for file in $(find ${alignmentDir}/filtered_bam -name "*.bam"); do
    sample=$(basename $file .bam)

    echo "Peak calling for sample ${sample}"
    macs2 callpeak -t ${alignmentDir}/filtered_bam/${sample}.bam \
        -f BAMPE \
        -g ${genomeSize} \
        -n ${sample} \
        -q ${qValue} \
        --nolambda \
        -B \
        --SPMR \
        --call-summits \
        --${cuttoff} \
        --outdir ${peakCallingDir}/macs2
done

# [END]