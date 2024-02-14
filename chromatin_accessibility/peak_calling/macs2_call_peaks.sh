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
echo "  Genome size: ${genomeSize}"
echo "  q-value: ${qValue}"
echo "  p-value: ${pValue}"
echo "  Perform cutoff analysis: ${cutoffAnalysis}"

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
case ${cutoffAnalysis} in
    true | True | TRUE | T | t )
        cutoffAnalysis="--cutoff-analysis"
        cutoffDir="/cutoff_analysis"
        if [ ! -d ${peakCallingDir}/macs2${cutoffDir} ]; then
            mkdir ${peakCallingDir}/macs2${cutoffDir}
        fi
        ;;
    *)
        cutoffAnalysis=""
        cutoffDir=""
        ;;
esac

# Determine if p-value is used to replace q-value for peak calling
# Determine is ${pValue} is set to a numeric value
if [[ ${pValue} =~ ^0\.[0-9]+$ ]]; then
    cutoff="-p ${pValue}"
    subDir="p${pValue}"
    echo "...Using p-value ${pValue} for peak calling"
    echo "...Overriding default q-value"
else
    cutoff="-q ${qValue}"
    subDir="q${qValue}"
fi

# Subdirectory for peak calling output based on parameters specified
subDir=${subDir}_nolambda

if [ ! -d ${peakCallingDir}/macs2${cutoffDir}/${subDir} ]; then
    mkdir ${peakCallingDir}/macs2${cutoffDir}/${subDir}
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
        ${cutoff} \
        --nolambda \
        -B \
        --SPMR \
        --call-summits \
        ${cutoffAnalysis} \
        --outdir ${peakCallingDir}/macs2${cutoffDir}/${subDir}
done

# [END]