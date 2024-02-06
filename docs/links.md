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
  ln ./transcriptome/pre-alignment_processing/fastqc_pre-trimming.sh \
    ./chromatin_accessibility/pre-alignment_processing/fastqc_pre-trimming.sh
  ```
  - `transcriptome` -> `protein_dna_interaction`
  
  ```bash
  ln ./transcriptome/pre-alignment_processing/fastqc_pre-trimming.sh \
    ./protein_dna_interaction/pre-alignment_processing/fastqc_pre-trimming.sh
  ```

- Trimming of reads with `fastp`

  - `chromatin_accessibility` -> `protein_dna_interaction`
  
  ```bash
  ln ./chromatin_accessibility/pre-alignment_processing/trimming_reads.sh \
    ./protein_dna_interaction/pre-alignment_processing/trimming_reads.sh
  ```

- FastQC processing of post-trimmed reads: 
  
  - `chromatin_accessibility` -> `protein_dna_interaction`
  
  ```bash
  ln ./chromatin_accessibility/pre-alignment_processing/fastqc_post-trimming.sh \
    ./protein_dna_interaction/pre-alignment_processing/fastqc_post-trimming.sh
  ```

- Alignment with Bowtie2
  
  - `transcriptome` -> `chromatin_accessibility`
  
  ```bash
  ln ./chromatin_accessibility/alignment/bowtie2_align.sh \
    ./protein_dna_interaction/alignment/bowtie2_align.sh
  ```

- Deduplication of aligned reads
  
  - `chromatin_accessibility` -> `protein_dna_interaction`
  
  ```bash
  ln ./chromatin_accessibility/post-alignment_processing/dedup_bam.sh \
    ./protein_dna_interaction/post-alignment_processing/dedup_bam.sh
  ```

- Fragment size distribution
    
  - `chromatin_accessibility` -> `protein_dna_interaction`
  
  ```bash
  ln ./chromatin_accessibility/post-alignment_processing/fragment_size.sh \
    ./protein_dna_interaction/post-alignment_processing/fragment_size.sh
  ```

  