# Plotting fragment size distribution, using the output generated from
# the 'fragment_size.sh' script.

# The inputs of this script must contain the following:
# 1. The fragment size distribution file, which is the output of the
#    'fragment_size.sh' script. This file contains two columns, the first
#    column is the fragment size, and the second column is the number of
#    fragments with that size.

# The output of this script is a set of PDF file, which contains the fragment
# size distribution plot for each mark.

# Usage: Rscript plot_frag_size.R <frag_size_file_dir> <plot_width> <plot_height>

# Check if the required arguments were provided
if (length(commandArgs(trailingOnly = TRUE)) != 1) {
    stop("Please provide the fragment size distribution file directory.")
}

# Check if the required R packages are installed, if not, install them
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
    install.packages("ggpubr")
}

# Define the global variables
fragSizeDir <- commandArgs(trailingOnly = TRUE)[1]
plotWidth <- commandArgs(trailingOnly = TRUE)[2]
plotWidth <- as.numeric(plotWidth)
plotHeight <- commandArgs(trailingOnly = TRUE)[3]
plotHeight <- as.numeric(plotHeight)

# Given the directory containing the fragment size distribution files,
# obtain the list of file names
fragSizeFiles <- list.files(fragSizeDir, pattern = "_fragment_size.txt$", full.names = TRUE)
samples <- gsub("_fragment_size.txt", "", basename(fragSizeFiles))

# Define a data frame to store the fragment size distribution for all samples
fragDist <- data.frame(Size = numeric(), Weight = numeric(), Sample = character())

# For each sample, read the fragment size distribution file and merge it to the data frame
for (sample in samples) {
    # Read the fragment size distribution file
    fragSizeFile <- paste0(fragSizeDir, "/", sample, "_fragment_size.txt")
    fragSize <- read.table(fragSizeFile, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    # Rename the columns
    colnames(fragSize) <- c("Size", "Weight")
    # Calculate weight
    fragSize$Weight <- fragSize$Weight / sum(fragSize$Weight)
    # Add the sample name column
    fragSize$Sample <- sample
    # Merge the fragment size distribution to the data frame
    fragDist <- rbind(fragDist, fragSize)
}

fragDist$Sample <- factor(fragDist$Sample, levels = samples)

# Plot the fragment size distribution as a violin plot, using ggplot2
library(ggplot2)
plot <- ggplot(fragDist, aes(x = Sample, y = Size, weight = Weight, fill = Sample)) +
        geom_violin(bw = 5) +
        scale_y_continuous(limits = c(0, maxLength), breaks = seq(0, maxLength, 50)) +
        scale_fill_brewer(palette = "Set1") +
        scale_color_brewer(palette = "Set1") +
        theme_bw() +
        ggpubr::rotate_x_text(angle = 45) +
        ylab("Fragment size (bp)") +
        xlab("")

# Save the plot as a PDF file
pdfFile <- paste0(outDir, "/", "frag_dist.pdf")
pdf(pdfFile, width = plotWidth, height = plotHeight)
print(plot)
dev.off()

# [END]