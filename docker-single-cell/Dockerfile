FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL maintainer="Jielin Yang"
LABEL email="jielin.yang@sickkids.ca"
LABEL version="1.0"
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
  cmake \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R --vanilla -e ' \
  install.packages(c("tidyverse", "dplyr", "devtools", "formatR", "remotes", "selectr", "caTools", "BiocManager", "BiocGenerics", "ggplot2", "rmarkdown", "httr", "knitr", "xaringan", "bookdown", "data.table", "tidyr", "shiny", "plotly", "kableExtra", "stringr", "stringi", "DT", "reshape2", "here", "pheatmap", "cowplot", "ggrepel", "readxl", "CVXR", "mixtools", "randomForest", "minerva", "igraph", "SCpubr", "reticulate", "sf", "scCustomize", "NMF", "Seurat", "SeuratObject", "Signac"), repos = "https://cran.rstudio.com", dependencies = TRUE); \
  setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/")); \
  install.packages(c("BPCells", "presto", "glmGamPoi")); \
  BiocManager::install(c("GEOquery", "DESeq2", "edgeR", "limma", "ComplexHeatmap", "biomaRt", "GEOmetadb", "BiocGenerics", "Biobase", "GSEABase", "GSVA", "BiocStyle", "fgsea", "BiocParallel", "GSA", "RCy3", "circlize", "chromVAR", "gprofiler2", "viper", "minet", "scater", "scran", "SingleCellExperiment", "DropletUtils", "SingleR", "celldex", "SpatialExperiment", "nnSVG", "DelayedArray", "DelayedMatrixStats", "lme4", "batchelor", "HDF5Array", "terra", "ggrastr", "harmony", "AUCell", "motifmatchr", "TFBSTools", "JASPAR2020", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86", "BSgenome.Mmusculus.UCSC.mm10", "EnsDb.Mmusculus.v79", "GenomicRanges", "GenomicFeatures", "AnnotationHub", "Biostrings", "ChIPseeker", "ATACseqQC")); \
  devtools::install_github("stjude/ChIPseqSpikeInFree"); \
  devtools::install_github("satijalab/seurat-data", quiet = TRUE); \
  devtools::install_github("satijalab/azimuth", quiet = TRUE); \
  devtools::install_github("satijalab/seurat-wrappers", quiet = TRUE); \
  devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories()); \
  ArchR::installExtraPackages(); \
  devtools::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE); \
  devtools::install_github("j-y26/bioToolkit"); \
  devtools::install_github("cole-trapnell-lab/monocle3"); \
  devtools::install_github("jinworks/CellChat") \
'

# Install Miniconda and set up Python environment
RUN R --vanilla -e ' \
  reticulate::install_miniconda(); \
  reticulate::conda_create("r-reticulate", python_version = "3.11"); \
  reticulate::conda_install("r-reticulate", c("numpy", "pandas", "scipy", "matplotlib", "seaborn", "scikit-learn", "statsmodels", "jupyter", "ipykernel", "tqdm", "pyarrow", "h5py", "scanpy", "squidpy", "anndata", "muon", "scvi-tools", "leidenalg", "louvain", "umap-learn", "harmonypy", "scrublet", "scanorama", "cellrank", "pynndescent", "spatialdata", "spatialde", "napari", "opencv-python-headless", "plotly", "bokeh", "networkx", "igraph", "adjustText", "tensorflow", "keras", "torch", "pytorch-lightning")); \
  try(system("~/.local/share/r-miniconda/envs/r-reticulate/bin/python -m ipykernel install --user --name r-reticulate --display-name \"Python (R-reticulate)\"", intern = TRUE), silent = TRUE) \
'

# Configure R to use the conda environment
RUN echo 'RETICULATE_PYTHON <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python"' >> /usr/local/lib/R/etc/Renviron.site && \
    echo 'RETICULATE_PYTHON_ENV <- "r-reticulate"' >> /usr/local/lib/R/etc/Renviron.site && \
    echo '# Configure reticulate to find Python\ntry({\n  library(reticulate)\n  use_condaenv("r-reticulate", required = TRUE)\n})' > /home/rstudio/.Rprofile

# Set working directory
WORKDIR /home/rstudio