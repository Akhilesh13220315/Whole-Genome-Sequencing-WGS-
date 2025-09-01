---
title: "End-to-End Whole Genome Sequencing (WGS) Pipeline"
author: "Bioinformatics Research Team"
date: "2025-09-01"
format:
  html:
    toc: true
    code-fold: true
    self-contained: true
  pdf:
    toc: true
    geometry: margin=1in
bibliography: references.bib
---

## Introduction

Whole Genome Sequencing (WGS) involves sequencing an organism's entire genome to identify genetic variants. This pipeline processes raw sequencing reads (FASTQ) to produce annotated variants (VCF) for downstream analysis. It follows GATK best practices for germline variant calling and is designed for reproducibility in a GitHub repository (`bioinformatics-pipelines`).

**Dataset**: We use open data from NCBI SRA (SRR030257, *Escherichia coli* paired-end reads, ~200MB) for demonstration due to its manageable size. For human WGS, adapt paths and reference to hg38.

**Tools**:
- **FastQC**: Quality control of raw reads.
- **Trimmomatic**: Adapter trimming and filtering.
- **BWA-MEM**: Alignment to reference genome.
- **SAMtools**: BAM processing.
- **Picard**: Duplicate marking.
- **GATK**: Variant calling and recalibration.
- **VEP**: Variant annotation.
- **MultiQC**: Aggregate QC reports.

**Setup**:
- Create a GitHub repository: `bioinformatics-pipelines`.
- Directory structure:
