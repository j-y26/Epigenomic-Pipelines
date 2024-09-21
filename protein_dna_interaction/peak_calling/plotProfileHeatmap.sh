# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Peak coverage matrix directory: ${peakCallingDir}/peak_coverage/${rawPeaks}"
echo "  TSS coverage matrix directory: ${peakCallingDir}/tss_coverage/${rawPeaks}"

# Iterate over the marks
for file in $(find ${markedSamples} -type f -name "samples_*.txt"); do
    mark=$(basename ${file} | sed 's/samples_//g' | sed 's/.txt//g')
    echo "Processing mark ${mark}"

    # Plotting profile plot and heatmap for peaks
    echo "Plotting profile plot and heatmap for peaks"
    plotProfile -m ${peakCallingDir}/peak_coverage/${rawPeaks}/${mark}_coverage_matrix_peaks.gz \
        -out ${peakCallingDir}/peak_coverage/${rawPeaks}/${mark}_profile_plot_peaks.pdf \
        --plotTitle "Average coverage over peaks" \
        --plotFileFormat "pdf" \
        --perGroup \
        --plotHeight 10 \
        --plotWidth 10 \
        --refPointLabel "Peak Center"

    plotHeatmap -m ${peakCallingDir}/peak_coverage/${rawPeaks}/${mark}_coverage_matrix_peaks.gz \
        -out ${peakCallingDir}/peak_coverage/${rawPeaks}/${mark}_heatmap_peaks.pdf \
        --plotTitle "Heatmap of coverage over peaks" \
        --plotFileFormat "pdf" \
        --perGroup \
        --refPointLabel "Peak Center"

    # Plotting profile plot and heatmap for TSS
    echo "Plotting profile plot and heatmap for TSS"
    plotProfile -m ${peakCallingDir}/tss_coverage/${rawPeaks}/${mark}_coverage_matrix_tss.gz \
        -out ${peakCallingDir}/tss_coverage/${rawPeaks}/${mark}_profile_plot_tss.pdf \
        --plotTitle "Average coverage over TSS" \
        --plotFileFormat "pdf" \
        --perGroup \
        --plotHeight 10 \
        --plotWidth 10 \
        --refPointLabel "TSS" \
        --startLabel "TSS - ${beforeTSSLength}" \
        --endLabel "TSS + ${beforeTSSLength}"

    plotHeatmap -m ${peakCallingDir}/tss_coverage/${rawPeaks}/${mark}_coverage_matrix_tss.gz \
        -out ${peakCallingDir}/tss_coverage/${rawPeaks}/${mark}_heatmap_tss.pdf \
        --plotTitle "Heatmap of coverage over TSS" \
        --plotFileFormat "pdf" \
        --perGroup \
        --refPointLabel "TSS" \
        --startLabel "TSS - ${beforeTSSLength}" \
        --endLabel "TSS + ${beforeTSSLength}"

    echo "Profile plot and heatmap generated for mark ${mark}"
done

# [END]