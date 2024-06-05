# Plotting the general statistics of the peak calling results, including
# the number of peaks and the peak width distribution.
# A bar plot is generated for all samples for the number of peaks.
# A density plot is generated for all samples for the peak width distribution.

# Usage: Rscript plot_peak_stats.R <peak_bed_dir> <output_dir> [<plot_width>] [<plot_height>]

# Check if the required arguments were provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_peak_stats.R <peak_bed_dir> <output_dir> [<plot_width>] [<plot_height>]")
}

# Define the global variables
peakBedDir <- commandArgs(trailingOnly = TRUE)[1]
if (! dir.exists(peakBedDir)) {
    stop(cat("The peak bed directory does not exist.\n",
             "Usage: Rscript plot_peak_stats.R <peak_bed_dir> <output_dir> [<plot_width>] [<plot_height>]"))
}
outputDir <- commandArgs(trailingOnly = TRUE)[2]
if (! dir.exists(outputDir)) {
    stop(cat("The output directory does not exist.\n",
             "Usage: Rscript plot_peak_stats.R <peak_bed_dir> <output_dir> [<plot_width>] [<plot_height>]"))
}

if (length(commandArgs(trailingOnly = TRUE)) == 4) {
    plotWidth <- commandArgs(trailingOnly = TRUE)[3]
    plotWidth <- as.numeric(plotWidth)
    plotHeight <- commandArgs(trailingOnly = TRUE)[4]
    plotHeight <- as.numeric(plotHeight)
} else {
    plotWidth <- 8
    plotHeight <- 6
}

# Check if the required R packages are installed, if not, install them
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
    install.packages("ggpubr")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
    install.packages("reshape2")
}
library(ggplot2)
library(ggpubr)
library(reshape2)

# Create the output directory if it does not exist
if (! dir.exists(outputDir)) {
    dir.create(outputDir)
}

# Summary of all the peak counts
peakCountSummary <- data.frame(Sample = character(), PeakCount = numeric())
maxPeakWidth90 <- 0

# Plot the fragment size distribution

# Get the sample names from all the peak bed files
peakBedFiles <- list.files(peakBedDir, pattern = ".bed$", full.names = TRUE)
samples <- gsub(".bed", "", basename(peakBedFiles))

# Summary of the peak distribution for all samples
summaryDist <- data.frame(Sample = character(),
                          Min = numeric(), Q1 = numeric(), Median = numeric(),
                          Mean = numeric(), Q3 = numeric(), Max = numeric())

# Get the peak count for each sample as well as the bed file
peakCounts <- c()
peakWidths <- list()
for (sample in samples) {
    peakBedFile <- paste0(peakBedDir, "/", sample, ".bed")
    if (! file.exists(peakBedFile)) {
        stop(cat("The peak bed file '", peakBedFile, "' does not exist.\n",
                    "Usage: Rscript plot_peak_stats.R <peak_bed_dir> <sample_matrix_csv> <output_dir> [<plot_width>] [<plot_height>]"))
    }

    # The number of peaks in the bed file
    peakCount <- as.numeric(system(paste("wc -l", peakBedFile, "| cut -d ' ' -f 1"), intern = TRUE))
    peakCounts <- c(peakCounts, peakCount)

    # The peak width distribution
    peakWidth <- read.table(peakBedFile, header = FALSE, stringsAsFactors = FALSE)
    peakWidth <- peakWidth$V3 - peakWidth$V2 + 1
    peakWidths[[sample]] <- peakWidth

    # Summary of the peak distribution
    summaryDist <- rbind(summaryDist, data.frame(Sample = sample,
                                                 Min = min(peakWidth), 
                                                 Q1 = quantile(peakWidth, 0.25),
                                                 Median = median(peakWidth), 
                                                 Mean = mean(peakWidth),
                                                 Q3 = quantile(peakWidth, 0.75), 
                                                 Max = max(peakWidth)))
    
    # Find the 90th percentile of the peak width, if max is too large
    if (max(peakWidth) > 10000) {
        peakWidth90 <- quantile(peakWidth, 0.90)
        maxPeakWidth90 <- max(maxPeakWidth90, peakWidth90)
    }
    
}

# Plotting the number of peaks for each sample
# Create a data frame for the peak counts
peakCountData <- data.frame(Sample = samples, PeakCount = peakCounts)
peakCountSummary <- rbind(peakCountSummary, peakCountData)

# Plot the peak count for each sample
peakCountPlot <- ggplot(peakCountData, aes(x = Sample, y = PeakCount)) +
    geom_bar(stat = "identity", fill = "#2c2494", width = 0.6, position = position_dodge(width = 0.8)) +
    geom_text(aes(label = PeakCount), vjust = -0.3, color = "black", size = 5) +
    labs(x = "Sample", y = "Number of Peaks") +
    theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
    theme_bw()
pdfFile <- paste0(outputDir, "/", "peak_count_bar.pdf")
pdf(pdfFile, width = plotWidth, height = plotHeight)
print(peakCountPlot)
dev.off()

# Plotting the peak width distribution for each sample
peakWidthData <- melt(peakWidths)
colnames(peakWidthData) <- c("PeakWidth", "Sample")
peakWidthData <- peakWidthData[peakWidthData$PeakWidth <= maxPeakWidth90, ]
cat("Plotting only peak widths less than ", maxPeakWidth90, ".\n")
peakWidthPlot <- ggplot(peakWidthData, aes(x = PeakWidth, color = Sample)) +
    geom_density() +
    labs(x = "Peak Width", y = "Density") +
    theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
    theme_bw()
pdfFile <- paste0(outputDir, "/", "peak_width_density.pdf")
pdf(pdfFile, width = plotWidth, height = plotHeight)
print(peakWidthPlot)
dev.off()

write.csv(summaryDist, paste0(outputDir, "/", "peak_width_summary.csv"), row.names = FALSE)


cat("Peak statistics plots have been generated successfully.\n")

write.csv(peakCountSummary, paste0(outputDir, "/peak_count_summary.csv"), row.names = FALSE)
cat("Peak statistics has been saved to '", paste0(outputDir, "/peak_count_summary.csv"), "'.\n")

# [END]