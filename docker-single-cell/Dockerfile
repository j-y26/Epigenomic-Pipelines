FROM bioconductor/bioconductor_docker:RELEASE_3_19
LABEL maintainer="Jielin Yang"
LABEL email="jielin.yang@sickkids.ca"
LABEL version="1.0"
LABEL description="Docker image for single-cell epigenomic analysis"

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
stringi stringr \
CVXR mixtools \
randomForest \
minerva \
igraph \
Seurat \
SCpubr \
Signac \
reticulate

# Install Bioconductor packages
RUN R --vanilla -e 'BiocManager::install(c("GEOQuery", \
"DESeq2", "edgeR", "limma", \
"ComplexHeatmap", "biomaRt","GEOmetadb","BiocGenerics", "Biobase"))'

RUN R --vanilla -e 'BiocManager::install(c("GSEABase", \
"GSVA", "BiocStyle", "fgsea", "BiocParallel"))'

RUN R --vanilla -e 'BiocManager::install(c("GSA", "RCy3", "seurat", \
"circlize", "chromVAR", "gprofiler2"))'

RUN R --vanilla -e 'BiocManager::install(c("viper", "minet"))'

# Install R packages from GitHub
RUN R --vanilla -e 'devtools::install_github("stjude/ChIPseqSpikeInFree")'

# Install miniconda
RUN R --vanilla -e 'reticulate::install_miniconda()'
