---
layout: default
title: Links
parent: Further Documentations
nav_order: 2
---

# Links used across pipelines

The use of hard links is to make sure that common tasks shared across
different pipelines are consistent and to reduce the risk of errors. The
following hard links are used across the pipelines:

Suppose `.` is the root of the project directory.

- FastQC processing of pre-trimmed reads: 
  
  - `transcriptome` -> `chromatin_accessibility`
  
  ```bash
  ln ./transcriptome/pre-alignment_processing/fastqc_pre-trimming.sh ./chromatin_accessibility/pre-alignment_processing/fastqc_pre-trimming.sh
  ```
  - `transcriptome` -> `protein_dna_interaction`
  
  ```bash
  ln ./transcriptome/pre-alignment_processing/fastqc_pre-trimming.sh ./protein_dna_interaction/pre-alignment_processing/fastqc_pre-trimming.sh
  ```

  


