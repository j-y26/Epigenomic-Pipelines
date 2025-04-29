FROM bioconductor/bioconductor_docker:RELEASE_3_21
LABEL maintainer="Jielin Yang"
LABEL email="jielin.yang@sickkids.ca"
LABEL version="1.0"
LABEL description="Docker image for single-cell and spatial epigenomic analysis"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    python3-dev \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN install2.r \
-d TRUE \
-r "https://cran.rstudio.com" \
RColorBrewer ggplot2 devtools \
rmarkdown httr knitr \
xaringan bookdown \
data.table \
tidyverse dplyr \
tidyr shiny \
plotly \
kableExtra \
stringr stringi \
DT \
reshape2 here pheatmap \
devtools cowplot \
ggrepel \
readxl \
CVXR mixtools \
randomForest \
minerva \
igraph \
SCpubr \
reticulate \
sf \
scCustomize \
NMF

# Single cell packages
RUN install2.r \
-d TRUE \
-r "https://cran.rstudio.com" \
Seurat Signac

# Add custom repositories for BPCells and other packages
RUN R --vanilla -e 'setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))'

# Install BPCells, presto, and glmGamPoi
RUN R --vanilla -e 'install.packages(c("BPCells", "presto", "glmGamPoi"))'

# Install Bioconductor packages for general genomics
RUN R --vanilla -e 'BiocManager::install(c("GEOQuery", \
  "DESeq2", "edgeR", "limma", \
  "ComplexHeatmap", "biomaRt","GEOmetadb","BiocGenerics", "Biobase"))'

# Install Bioconductor packages for pathway analysis
RUN R --vanilla -e 'BiocManager::install(c("GSEABase", \
  "GSVA", "BiocStyle", "fgsea", "BiocParallel"))'

# Install additional Bioconductor packages
RUN R --vanilla -e 'BiocManager::install(c("GSA", "RCy3", \
  "circlize", "chromVAR", "gprofiler2"))'

# Install network analysis packages
RUN R --vanilla -e 'BiocManager::install(c("viper", "minet"))'

# Install additional single-cell and spatial packages from Bioconductor
RUN R --vanilla -e 'BiocManager::install(c("scater", "scran", "SingleCellExperiment", \
  "DropletUtils", "SingleR", "celldex", "SpatialExperiment", \
  "nnSVG", "DelayedArray", "DelayedMatrixStats", "lme4", "batchelor", "HDF5Array", "terra", "ggrastr", \
  "harmony", "AUCell"))'

# Install packages for chromatin accessibility analysis
RUN R --vanilla -e 'BiocManager::install(c("motifmatchr", "TFBSTools", "JASPAR2020", \
  "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86", \
  "ChIPseeker", "ATACseqQC", "ChromVAR"))'

# Install R packages from GitHub
RUN R --vanilla -e 'devtools::install_github("stjude/ChIPseqSpikeInFree")'
RUN R --vanilla -e 'devtools::install_github("satijalab/seurat-data", quiet = TRUE)'
RUN R --vanilla -e 'devtools::install_github("satijalab/azimuth", quiet = TRUE)'
RUN R --vanilla -e 'devtools::install_github("satijalab/seurat-wrappers", quiet = TRUE)'
RUN R --vanilla -e 'devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())'
RUN R --vanilla -e 'ArchR::installExtraPackages()'

# Install spatial packages
RUN R --vanilla -e 'devtools::install_github("cole-trapnell-lab/monocle3")'
RUN R --vanilla -e 'devtools::install_github("jinworks/CellChat")'

# Install miniconda for Python integration
RUN R --vanilla -e 'reticulate::install_miniconda()'

# Create Python environment with specific version
RUN R --vanilla -e 'reticulate::conda_create("r-reticulate", python_version = "3.10")'

# Install core Python packages for data science and single-cell analysis
RUN R --vanilla -e 'reticulate::conda_install("r-reticulate", c("numpy", "pandas", "scipy", "matplotlib", "seaborn", "scikit-learn", "statsmodels", "jupyter", "ipykernel", "tqdm", "pyarrow", "h5py"))'

# Install single-cell specific Python packages
RUN R --vanilla -e 'reticulate::conda_install("r-reticulate", c("scanpy", "squidpy", "anndata", "muon", "scvi-tools", "leidenalg", "louvain", "umap-learn", "harmonypy", "scrublet", "scanorama", "cellrank", "pynndescent"))'

# Install spatial analysis packages
RUN R --vanilla -e 'reticulate::conda_install("r-reticulate", c("squidpy", "spatialdata", "spatialde", "napari", "opencv-python-headless"))'

# Install plotting and visualization packages
RUN R --vanilla -e 'reticulate::conda_install("r-reticulate", c("plotly", "bokeh", "networkx", "igraph", "adjustText"))'

# Install machine learning packages
RUN R --vanilla -e 'reticulate::conda_install("r-reticulate", c("tensorflow", "keras", "torch", "pytorch-lightning", "scvi-tools"))'

# Make Python environment accessible via Jupyter
RUN R --vanilla -e 'reticulate::conda_install("r-reticulate", c("jupyter", "ipykernel"))'
RUN R --vanilla -e 'system("~/.local/share/r-miniconda/envs/r-reticulate/bin/python -m ipykernel install --user --name r-reticulate --display-name \"Python (R-reticulate)\"")' || true

# Configure R to use the conda environment by default
RUN echo 'RETICULATE_PYTHON <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python"' >> /usr/local/lib/R/etc/Renviron.site && \
    echo 'RETICULATE_PYTHON_ENV <- "r-reticulate"' >> /usr/local/lib/R/etc/Renviron.site

# Create R initialization script to verify Python integration
RUN echo '# Configure reticulate to find Python\ntry({\n  library(reticulate)\n  use_condaenv("r-reticulate", required = TRUE)\n})' > /home/rstudio/.Rprofile

# Set working directory
WORKDIR /home/rstudio