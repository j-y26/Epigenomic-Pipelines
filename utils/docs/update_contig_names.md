---
layout: default
title: Chromosome Name Conversion
parent: Utilities
nav_order: 6
---

# Converting Chromosome Names

Here we make use of the useful resources from the Github repository
[ChromosomeMappings](https://github.com/dpryan79/ChromosomeMappings) to convert
chromosome names between different genome assemblies.

## Installation

```bash
conda install -c bioconda cvbio
```

To update any delimited files (e.g. BED, GTF, VCF) with the new chromosome names,
use the following command:

```bash
wget https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/<mapping_file>

cvbio UpdateContigNames \
    -i <input_file> \
    -o <output_file> \
    -m <mapping_file> \
    --comment-char '#' \
    --columns 0 \
    --skip-missing false
```