---
layout: default
title: Parsing GTF
parent: Utilities
nav_order: 2
---

# Parsing GTF

We have provided several utility scripts to parse GTF files into tab-delimited 
files, which allows users to easily use table operations (e.g., working as dataframes
in R or pandas in Python) to manipulate the data.

The scripts will extract either the gene, transcript, or exon information from 
the GTF file.

- Gene information: [`extract_gene_from_gtf.py`](https://github.com/j-y26/Epigenomic-Pipelines/blob/535376041012783fb84b3b3c1e58f0811c13877e/utils/extract_gene_from_gtf.py)
- Transcript information: [`extract_transcript_from_gtf.py`](https://github.com/j-y26/Epigenomic-Pipelines/blob/535376041012783fb84b3b3c1e58f0811c13877e/utils/extract_transcript_from_gtf.py)
- Exon information: [`extract_exon_from_gtf.py`](https://github.com/j-y26/Epigenomic-Pipelines/blob/535376041012783fb84b3b3c1e58f0811c13877e/utils/extract_exon_from_gtf.py)

To make use of the scripts, users can run the following command:

```bash
python3 extract_gene_from_gtf.py <path_to_gtf_file> <output_file>
```

The output file will be a tab-delimited file containing gene/transcript/exon 
attributes.