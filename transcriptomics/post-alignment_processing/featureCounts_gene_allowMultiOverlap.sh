# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "nameSortedDir: ${nameSortedDir}"
echo "featureCountsDir: ${featureCountsDir}"
echo "outCountFile: ${outCountFile}"

# [Main]

# Create output directory if it doesn't exist
if [ ! -d ${featureCountsDir} ]; then
    mkdir ${featureCountsDir}
fi

# Count the number of reads that map to each gene using featureCounts
# Here, we generate the count by gene, for gene-based differential analysis

samples=$(find ${nameSortedDir} -type f -name '*_name_sorted.bam')
echo "Counting reads ..."
featureCounts \
  -T ${threads} \
  -p --countReadPairs \
  -t exon \
  -g gene_id \
  -O \
  -a ${genomeGtfFile} \
  -o ${featureCountsDir}/${outCountFile} \
  ${samples}
echo "Counting completed"

# [END] 