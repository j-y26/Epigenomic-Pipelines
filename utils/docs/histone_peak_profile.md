---
layout: default
title: Histone Peak Profiles for ChIP
parent: Utilities
nav_order: 5
---

# Histone Peak Profiles for ChIP

Histone marks have unique distributions across the genome, and are often have
different functions in gene regulation. Therefore, their distribution profiles
affects how peak callers should be used to identify peaks.

Here we summarize the distribution profiles of some common histone marks, many of which are adapted from the [ENCODE project](https://www.encodeproject.org/chip-seq/histone/).

In addition, RNA polymerase often has different post-translational
modifications on its tail to regulate RNA polymerase activity. These
modifications can be used to identify active transcription sites.
Phosphorylation of RNA polymerase II at serine 5 (Pol2S5P) is a common mark
for active transcription initiation, while phosphorylation at serine 2
(Pol2S2P) is a common mark for active transcription elongation. We also
include structural proteins that has stable binding sites across the genome.

Here we summarize the distribution profiles of some common histone marks and RNA polymerase marks:

| Broad peak | Narrow peak | Exception |
|------------|-------------|-----------|
H3F3A | H2AFZ | H3K9me3
H3K27me3 | H3ac | 
H3K36me3 | H3K27ac |
H3K4me1 | H3K4me2 |
H3K79me2 | H3K4me3 |  
H3K79me3 | H3K9ac |
H3K9me1 | H2AUb119 |
H3K9me2 | CTCF |
H4K20me1 | RNAP-pSer5 |
RNAP-pSer2 | |

Although H3K9me3 is an exception due to the fact that most of these peaks reside in repetitive regions, it can be called as a broad peak.