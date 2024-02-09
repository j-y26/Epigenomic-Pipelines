# Plotting fragment size distribution, using the output generated from
# the 'fragment_size.sh' script.

# The inputs of this script must contain the following:
# 1. The fragment size distribution file, which is the output of the
#    'fragment_size.sh' script. This file contains two columns, the first
#    column is the fragment size, and the second column is the number of
#    fragments with that size.
# 2. A sample matrix file, containing a 'Mark' column and a 'Label' column.
#    The 'Mark' column contains the name of the target, and the 'Label'
#    column contains the label of the sample, which corresponds to the file
#    name of the fragment size distribution file.

# The output of this script is a set of PDF file, which contains the fragment 
# size distribution plot for each mark.

# Usage: Rscript plot_frag_size.R <frag_size_file_dir> <sample_matrix> <max_fragment_length> [<plot_width>] [<plot_height>]

# Check if the required arguments were provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_frag_size.R <frag_size_file_dir> <sample_matrix> <max_fragment_length> [<plot_width>] [<plot_height>]")
}

# Define the global variables
fragSizeDir <- commandArgs(trailingOnly = TRUE)[1]
if (! dir.exists(fragSizeDir)) {
    stop(cat("The fragment size distribution directory does not exist.\n",
             "Usage: Rscript plot_frag_size.R <frag_size_file_dir> <sample_matrix> <max_fragment_length> [<plot_width>] [<plot_height>]"))
}
sampleMatrix <- commandArgs(trailingOnly = TRUE)[2]
if (! file.exists(sampleMatrix)) {
    stop(cat("The sample matrix file does not exist.\n",
             "Usage: Rscript plot_frag_size.R <frag_size_file_dir> <sample_matrix> <max_fragment_length> [<plot_width>] [<plot_height>]"))
}
maxLength <- commandArgs(trailingOnly = TRUE)[3]
maxLength <- as.numeric(maxLength)

if (length(commandArgs(trailingOnly = TRUE)) == 5 {
    plotWidth <- commandArgs(trailingOnly = TRUE)[4]
    plotWidth <- as.numeric(plotWidth)
    plotHeight <- commandArgs(trailingOnly = TRUE)[5]
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

# Read the sample matrix file
sampleMatrix <- read.csv(sampleMatrix, header = TRUE, stringsAsFactors = FALSE)

# Check if sample matrix contains the required columns
if (!all(c("Mark", "Label") %in% colnames(sampleMatrix))) {
    stop("The sample matrix file must contain 'Mark' and 'Label' columns, where 'Mark' contains the name of the target, and 'Label' contains the label of the sample.")
}

# Get the list of marks
marks <- unique(sampleMatrix$Mark)

# Plot the fragment size distribution for each mark
for (mark in marks) {
    # Get the sample names for the current mark
    samples <- sampleMatrix$Label[sampleMatrix$Mark == mark]

    # Define a data frame to store the fragment size distribution for the current mark
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
            scale_y_continuous(limits = c(0, 750), breaks = seq(0, maxLength, 50)) +
            scale_fill_brewer(palette = "Set1") +
            scale_color_brewer(palette = "Set1") +
            theme_bw() +
            ggpubr::rotate_x_text(angle = 45) +
            ylab("Fragment size (bp)") +
            xlab("")

    # Save the plot as a PDF file
    pdfFile <- paste0(outDir, "/", mark, "_frag_dist.pdf")
    pdf(pdfFile, width = 8, height = 6)
    print(plot)
    dev.off()
}

# [END]