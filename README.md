[![](/docs/assets/page_logo.png)](https://j-y26.github.io/Epigenomic-Pipelines/)

### Table of Contents

- [Introduction](#introduction)
- [Implemented Pipelines](#implemented-pipelines)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

This repository contains a collection of pipelines for processing epigenomic NGS
data, including the analysis of the transcriptome, chromatin accessibility, and
protein-DNA interactions.

The pipelines primarily focus on the processing of raw sequencing data to the
generation of high-quality data files that can be used for downstream analysis.
The pipelines are designed to be modular, allowing for the easy addition of new
tools and methods.

## Implemented Pipelines

Currently, the following analysis pipelines are implemented:

- Transcriptome analysis (RNA-seq)
- Chromatin accessibility analysis (ATAC-seq)
- Protein-DNA interaction analysis (ChIP-seq, CUT&RUN, CUT&Tag)

Each pipeline is implemented as a set of modular steps, under a common
framework. The pipelines are implemented primarily using shell scripts, where
each step is implemented as a separate script. Ensure that you are using a
compatible shell (e.g. `bash`) and are using a Unix-based operating system.

## Installation

### Linux

To install the pipelines, clone the repository and install the required
dependencies.

```bash
git clone https://github.com/j-y26/Epigenomic-Pipelines.git
cd Epigenomic-Pipelines
```

In this pipeline, we assume that the required tools are installed and available
in the system path. The required tools are listed in the `requirements.md` file,
and all tools are publicly available.

### Docker - Linux for Bulk Epigenomic Processing

Alternatively, a Docker image is provided that contains all the required tools
and dependencies. To use the Docker image, ensure that Docker is installed and
run the following command:

```bash
docker run -it --user $(id -u):$(id -g) -v ${pwd}:/home jyang26/epigenomic-pipelines:latest
```

This will download the Docker image and start an interactive shell session.
Ensure you are in the correct directory where you want to mount the data or 
change the `${pwd}` to the correct directory.

To return to the same container once the container is stopped, run the following
command:

```bash
# Replace <container_id> with the container ID of the corresponding container
docker start -i <container_id>
docker attach <container_id>
```

### Docker - Single-cell Epigenomic Analysis

The docker image for single-cell analysis is based on Bioconductor releases.
Therefore, the image will build a RStudio server environment with version-controlled
packages installed. To pull and install the container:

```bash
docker run -e PASSWORD=changeit -v ${pwd}:/home/rstudio/projects -p 8787:8787 jyang26/single-cell-epigenomics:v1.0
```

A Python environment (Miniconda) is also installed in this environment. CLI of this
container can be accessed in a similar way.

## Usage

To use the pipelines, detailed descriptions and instructions are provided:

- [Transcriptome Analysis](transcriptome/README.md)
- [Chromatin Accessibility Analysis](chromatin_accessibility/README.md)
- [Protein-DNA Interaction Analysis](protein_dna_interaction/README.md)

Each directory contains a `README.md` file with detailed instructions on how to
use the pipeline. These directories also contain the scripts and configuration
files required to run the pipeline.

For each pipeline, we have provided an example configuration file that is
necessary to define the variables and parameters for the analysis. Importantly,
users must specify the path to the data input files in the configuration file,
since the bash scripts themselves are designed to be modular and generalizable
to any dataset.

Users can move the configuration files to any directory and make modifications
tailored to their specific data. However, ensure you have the correct path to
the configuration file when running the pipeline. Ensure you use the bash
syntax when working with the configuration file.

A general example for executing the scripts is as follows:

```bash
./epigenomic-pipelines/pipeline_folder/script.sh /path/to/config.sh
```

To reduce the risk of errors, some scripts that perform common tasks are created
using symbolic links. These scripts can be used the same way as the original
scripts. Documentation of these symbolic links can be found [here](docs/symlinks.md).

## Contributing

Contributions are welcome! For major changes, please open an issue first to
discuss what you would like to change. For any issues, please open an issue in
the [issue tracker](https://github.com/j-y26/Epigenomic-Pipelines/issues).

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE)
file for details.
