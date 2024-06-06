# Plotting the general statistics of the peak calling results
# A bar plot is generated for all samples with the same mark for the number of peaks.

# Usage: Rscript plot_peak_count_consensus.R <replicate_peak_dir> <consensus_peak_dir> <sample_matrix_csv> [<plot_width>] [<plot_height>]

# Check if the required arguments were provided
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) < 3) {
    stop("Usage: Rscript plot_peak_count_consensus.R <replicate_peak_dir> <consensus_peak_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]")
}

# Define the global variables
replicatePeakDir <- commandArgs(trailingOnly = TRUE)[1]
if (! dir.exists(replicatePeakDir)) {
    stop(cat("The replicate peak directory does not exist.\n",
             "Usage: Rscript plot_peak_count_consensus.R <replicate_peak_dir> <consensus_peak_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]"))
}
consensusPeakDir <- commandArgs(trailingOnly = TRUE)[2]
if (! dir.exists(consensusPeakDir)) {
    stop(cat("The consensus peak directory does not exist.\n",
             "Usage: Rscript plot_peak_count_consensus.R <replicate_peak_dir> <consensus_peak_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]"))
}
sampleMatrix <- commandArgs(trailingOnly = TRUE)[3]
if (! file.exists(sampleMatrix)) {
    stop(cat("The sample matrix file does not exist.\n",
             "Usage: Rscript plot_peak_count_consensus.R <replicate_peak_dir> <consensus_peak_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]"))
}

if (length(commandArgs(trailingOnly = TRUE)) == 5) {
    plotWidth <- commandArgs(trailingOnly = TRUE)[4]
    plotWidth <- as.numeric(plotWidth)
    plotHeight <- commandArgs(trailingOnly = TRUE)[5]
    plotHeight <- as.numeric(plotHeight)
} else {
    plotWidth <- 10
    plotHeight <- 6
}

# Check if the required R packages are installed, if not, install them
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
    install.packages("ggpubr")
}
library(ggplot2)
library(ggpubr)

# Read the sample matrix file
sampleMatrix <- read.csv(sampleMatrix, header = TRUE, stringsAsFactors = FALSE)
sampleMatrix$Peak_Type <- tolower(sampleMatrix$Peak_Type)
sampleMatrix <- sampleMatrix[sampleMatrix$Peak_Type %in% c("narrow", "broad"), ]

# Check if sample matrix contains the required columns
if (! all(c("Mark", "Label") %in% colnames(sampleMatrix))) {
    stop("The sample matrix file must contain 'Mark' and 'Label' columns.")
}

# Get the list of marks
marks <- unique(sampleMatrix$Mark)

# Plot the fragment size distribution for each mark
for (mark in marks) {
    # Get the sample names for the current mark
    samples <- sampleMatrix$Label[sampleMatrix$Mark == mark]
    samples <- samples[order(samples)]
    
    # Get the peak count for each sample as well as the bed file
    peakCounts <- c()
    for (sample in samples) {
        replicatePeakFile <- paste0(replicatePeakDir, "/", sample, ".bed")
        if (! file.exists(replicatePeakFile)) {
            stop(cat("The replicate peak file '", replicatePeakFile, "' does not exist.\n",
                     "Usage: Rscript plot_peak_count_consensus.R <replicate_peak_dir> <consensus_peak_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]"))
        }
        replicatePeakCount <- as.numeric(system(paste("wc -l", replicatePeakFile, "| cut -d ' ' -f 1"), intern = TRUE))
        peakCounts <- c(peakCounts, replicatePeakCount)
    }

    # Find the consensus peak file, which could be multiple bed files containing the same mark in the file name
    consensusPeakFilesCov <- list.files(consensusPeakDir, pattern = paste0(mark, "_consensus_cov.bed"), full.names = TRUE)
    consensusPeakFilesOverlap <- list.files(consensusPeakDir, pattern = paste0(mark, "_consensus_overlap.bed"), full.names = TRUE)

    # The number of peaks in the consensus peak file
    consensusSamples <- c()
    for (consensusPeakFile in consensusPeakFilesCov) {
        peakCounts <- c(peakCounts, as.numeric(system(paste("wc -l", consensusPeakFile, "| cut -d ' ' -f 1"), intern = TRUE)))
        consensusSamples <- c(consensusSamples, basename(consensusPeakFile))
    }
    for (consensusPeakFile in consensusPeakFilesOverlap) {
        peakCounts <- c(peakCounts, as.numeric(system(paste("wc -l", consensusPeakFile, "| cut -d ' ' -f 1"), intern = TRUE)))
        consensusSamples <- c(consensusSamples, basename(consensusPeakFile))
    }
    
    # Create a data frame for the peak counts
    peakCountData <- data.frame(Sample = c(samples, consensusSamples), PeakCount = peakCounts)
    
    # Plot the peak count for each sample
    peakCountPlot <- ggplot(peakCountData, aes(x = Sample, y = PeakCount, fill = Sample)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        geom_text(aes(label = PeakCount), vjust = -0.3, color = "black", size = 5) +
        labs(x = "Sample", y = "Number of Peaks") +
        theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16)) +
        theme_bw()
    pdfFile <- paste0(consensusPeakDir, "/", mark, "_peak_count_bar.pdf")
    pdf(pdfFile, width = plotWidth, height = plotHeight)
    print(peakCountPlot)
    dev.off()

    print(paste("The peak count bar plot for", mark, "was saved to", pdfFile))
}


# [END]