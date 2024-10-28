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

We provided a script to parse any tab-delimited file with the
chromosome names. The script is located at [convert_contig_names.py](../convert_contig_names.py)

```bash
python3 convert_contig_names.py <input_file> <mapping_file> <output_file> [<chromosome_column>]
```

By default, if the column number is not specified, the script will assume the
first column contains the chromosome names. Alternatively, the column number
can be specified as the fourth argument, which should be 1-indexed, i.e., the
first column is 1.