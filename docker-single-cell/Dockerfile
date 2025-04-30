FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL maintainer="Jielin Yang"
LABEL email="jielin.yang@sickkids.ca"
LABEL version="1.1"
LABEL description="Docker image for single-cell and spatial epigenomic analysis"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
  libxml2-dev \
  libcairo2-dev \
  libsqlite3-dev \
  libmariadbd-dev \
  libmariadb-dev \
  libpq-dev \
  libssh2-1-dev \
  unixodbc-dev \
  libsasl2-dev \
  libhdf5-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libpng-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  python3-dev \
  python3-pip \
  libbz2-dev \
  zlib1g-dev \
  gfortran \
  libblas-dev \
  liblapack-dev \
  libgl1-mesa-dev \
  libglu1-mesa-dev \
  texlive \
	texlive-latex-extra \
	texlive-fonts-extra \
	texlive-bibtex-extra \
	texlive-science \
	texi2html \
	texinfo \
  cmake \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/miniconda && \
    rm /tmp/miniconda.sh && \
    source /opt/miniconda/bin/activate

# Install R packages
RUN R --vanilla -e ' \
  options(repos = c(CRAN = "https://cran.rstudio.com")); \
  install.packages(c("tidyverse", "dplyr", "devtools", "formatR", "remotes", "selectr", "caTools", "ggplot2", "rmarkdown", "httr", "knitr", "xaringan", "bookdown", "data.table", "tidyr", "shiny", "plotly", "kableExtra", "stringr", "stringi", "DT", "reshape2", "here", "pheatmap", "cowplot", "ggrepel", "readxl", "writexl", "readr", "CVXR", "mixtools", "randomForest", "minerva", "igraph", "ggfortify", "rlang", "sf", "reticulate"), repos = "https://cran.rstudio.com", dependencies = TRUE) \
'

RUN R --vanilla -e ' \
  options(repos = c(CRAN = "https://cran.rstudio.com")); \
  BiocManager::install(c("AnnotationDbi", "AUCell", "clusterProfiler", "ComplexHeatmap", "decoupleR", "enrichplot", "infercnv", "Nebulosa", "org.Hs.eg.db", "org.Mm.eg.db", "UCell", "BiocFileCache", "dittoSeq", "scuttle", "DelayedArray", "IRanges", "MAST", "monocle", "rtracklayer", "S4Vectors", "SummarizedExperiment", "Rsamtools", "BSgenome", "biovizBase", "Biostrings", "multtest", "glmGamPoi", "GEOquery", "DESeq2", "edgeR", "limma", "biomaRt", "GEOmetadb", "BiocGenerics", "Biobase", "GSEABase", "GSVA", "BiocStyle", "fgsea", "BiocParallel", "GSA", "RCy3", "circlize", "chromVAR", "gprofiler2", "viper", "minet", "scater", "scran", "SingleCellExperiment", "DropletUtils", "SingleR", "celldex", "SpatialExperiment", "nnSVG", "DelayedMatrixStats", "lme4", "batchelor", "HDF5Array", "terra", "ggrastr", "harmony", "motifmatchr", "TFBSTools", "JASPAR2020", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86", "BSgenome.Mmusculus.UCSC.mm10", "EnsDb.Mmusculus.v79", "GenomicRanges", "GenomicFeatures", "AnnotationHub", "Biostrings", "ChIPseeker", "ATACseqQC"), dependencies = TRUE); \
  devtools::install_github("stjude/ChIPseqSpikeInFree"); \
  devtools::install_github("satijalab/seurat-data"); \
  devtools::install_github("satijalab/azimuth"); \
  devtools::install_github("satijalab/seurat-wrappers"); \
  devtools::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE); \
  devtools::install_github("j-y26/bioToolkit"); \
  devtools::install_github("cole-trapnell-lab/monocle3"); \
  devtools::install_github("immunogenomics/presto"); \
  devtools::install_github("bnprks/BPCells/r"); \
  devtools::install_github("jinworks/CellChat") \
'

RUN R --vanilla -e ' \
  options(repos = c(CRAN = "https://cran.rstudio.com")); \
  install.packages(c("SCpubr", "scCustomize", "NMF", "Seurat", "SeuratObject", "Signac", "enrichR"), repos = "https://cran.rstudio.com", dependencies = TRUE) \
'