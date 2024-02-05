---
layout: default
title: Resources
parent: Utilities
nav_order: 1
---

# Resources

## Publicly Available Data

Some useful data are provided as resources. These resources are kept in its 
original format and are not modified. Please refer to the original publishing 
source for more information about the data.

- [mm10.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes)
- [hg38.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)
- [rn7.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.chrom.sizes)
- [mm10-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)
- [hg38-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)

## Compiled Resources

For convenience, some resources are compiled and provided here.

For interactively querying the biomaRt database in R, how different values are
encoded in the database is presented here:

- [Interactive biomaRt database index in R](./BioMart_Query.html)

Normalizing reads when generating bigWig files is a common practice. An
effective genome size is sometimes required. Depending on whether multi-mapping
reads are included, the effective genome size can be calculated differently.
DeepTools has provided some values to commonly used genomes. This information
can be found here:

- [Effective Genome Sizes](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html)

When **multi-mapping reads are included**, we use the non-N bases in the genome:

| Genome | Effective Genome Size |
|----------|----------|
| mm10/GRCm38 | `2652783500` |
| hg38/GRCh38 | `2913022398` |

When **multi-mapping reads are excluded**, we use the uniquely mappable genome size:

| Read Length | mm10/GRCm38 | hg38/GRCh38 |
|-------------|--------------|-------------|
| 50 | `2308125299` | `2701495711` |
| 75 | `2407883243` | `2747877702` |
| 100 | `2467481008` | `2805636231` |
| 150 | `2494787038` | `2862010428` |
| 200 | `2520868989` | `2887553103` |
| 250 | `2538590322` | `2898802627` |