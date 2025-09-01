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

bioinformatics-pipelines/ ├── data/raw/ ├── ref/ ├── results/ ├── scripts/ ├── README.md └── wgs-pipeline.qmd

- Install Quarto for rendering: `quarto install`.
- Install dependencies or use Docker (e.g., `nf-core/sarek` image).
- Render with: `quarto render wgs-pipeline.qmd`.

**Assumptions**:
- Input: Paired-end FASTQ files.
- Reference: *E. coli* genome (GCF_000005845.2).
- Environment: Linux/Colab with sufficient memory (~16GB).

## Step 1: Raw Read Quality Control

Quality control assesses raw FASTQ files for sequencing errors, adapter content, and base quality.

### Code
```{bash}
# Create output directory
mkdir -p results/fastqc

# Run FastQC
fastqc -o results/fastqc data/raw/SRR030257_1.fastq.gz data/raw/SRR030257_2.fastq.gz

Explanation





FastQC generates HTML reports for metrics like per-base quality, GC content, and adapter contamination.



Check for low-quality bases (Phred < 20) or overrepresented sequences.



Output: results/fastqc/*.html. Review for issues before proceeding.

Step 2: Adapter Trimming and Preprocessing

Trim adapters and low-quality bases to improve alignment accuracy.

Code

# Create output directory
mkdir -p results/trimmed

# Run Trimmomatic for paired-end reads
trimmomatic PE -phred33 \
  data/raw/SRR030257_1.fastq.gz data/raw/SRR030257_2.fastq.gz \
  results/trimmed/SRR030257_1_trimmed.fastq.gz results/trimmed/SRR030257_1_unpaired.fastq.gz \
  results/trimmed/SRR030257_2_trimmed.fastq.gz results/trimmed/SRR030257_2_unpaired.fastq.gz \
  ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Explanation





Trimmomatic removes adapters (e.g., TruSeq3) and low-quality bases.



Parameters:





ILLUMINACLIP: Removes adapter sequences.



LEADING/TRAILING: Trims bases with Phred < 3.



SLIDINGWINDOW: Removes windows with average Phred < 15.



MINLEN: Discards reads < 36bp.



Output: Paired (_trimmed) and unpaired reads.



Re-run FastQC on trimmed files to verify improvements.

Step 3: Download Reference and Data

Download open data and reference genome for E. coli.

Code

# Create directories
mkdir -p data/raw ref

# Download SRA data
fasterq-dump SRR030257 --split-files -O data/raw

# Download E. coli reference
wget -O ref/ecoli.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip ref/ecoli.fa.gz

# Index reference
bwa index ref/ecoli.fa
samtools faidx ref/ecoli.fa

Explanation





fasterq-dump: Downloads SRR030257 from NCBI SRA.



Reference: E. coli genome (GCF_000005845.2).



Indexing: Prepares reference for alignment (BWA) and BAM processing (SAMtools).

Step 4: Alignment to Reference Genome

Align trimmed reads to the reference genome using BWA-MEM.

Code

# Create output directory
mkdir -p results/aligned

# Align with BWA-MEM
bwa mem -R "@RG\tID:SRR030257\tSM:SRR030257\tLB:lib1\tPL:illumina" \
  ref/ecoli.fa results/trimmed/SRR030257_1_trimmed.fastq.gz results/trimmed/SRR030257_2_trimmed.fastq.gz \
  > results/aligned/SRR030257.sam

# Convert to BAM, sort, and index
samtools view -bS results/aligned/SRR030257.sam | samtools sort -o results/aligned/SRR030257_sorted.bam
samtools index results/aligned/SRR030257_sorted.bam

Explanation





BWA-MEM: Aligns paired-end reads to the reference.



Read Group (-R): Adds metadata for downstream tools.



SAMtools:





Converts SAM to BAM (view -bS).



Sorts BAM by coordinate (sort).



Indexes BAM for fast access (index).



Output: SRR030257_sorted.bam and .bai index.

Step 5: Post-Alignment Processing

Clean BAM files by marking duplicates and recalibrating base quality scores.

Code

# Mark duplicates with Picard
java -jar /path/to/picard.jar MarkDuplicates \
  I=results/aligned/SRR030257_sorted.bam \
  O=results/aligned/SRR030257_dedup.bam \
  M=results/aligned/SRR030257_metrics.txt

# Index deduplicated BAM
samtools index results/aligned/SRR030257_dedup.bam

# Base Quality Score Recalibration (BQSR) with GATK
gatk BaseRecalibrator \
  -I results/aligned/SRR030257_dedup.bam \
  -R ref/ecoli.fa \
  --known-sites ref/known.vcf \
  -O results/aligned/recal_data.table

gatk ApplyBQSR \
  -I results/aligned/SRR030257_dedup.bam \
  -R ref/ecoli.fa \
  --bqsr-recal-file results/aligned/recal_data.table \
  -O results/aligned/SRR030257_recal.bam

Explanation





Picard MarkDuplicates: Identifies and marks PCR duplicates to prevent bias in variant calling.



GATK BaseRecalibrator: Builds a recalibration model based on known variants (assumes known.vcf exists; for E. coli, may skip or use a public database like dbSNP).



ApplyBQSR: Adjusts base quality scores to improve variant calling accuracy.



Output: SRR030257_recal.bam (recalibrated BAM).

Step 6: Variant Calling

Identify germline variants using GATK HaplotypeCaller.

Code

# Create output directory
mkdir -p results/variants

# Run HaplotypeCaller
gatk HaplotypeCaller \
  -R ref/ecoli.fa \
  -I results/aligned/SRR030257_recal.bam \
  -O results/variants/SRR030257.g.vcf \
  -ERC GVCF

Explanation





HaplotypeCaller: Calls SNPs and indels in GVCF mode for scalability.



-ERC GVCF: Produces a genomic VCF for joint genotyping if multiple samples are used.



Output: SRR030257.g.vcf.

For multi-sample analysis, use GenotypeGVCFs (not shown for single sample).

Step 7: Variant Filtering

Filter variants to remove low-quality calls.

Code

gatk VariantFiltration \
  -R ref/ecoli.fa \
  -V results/variants/SRR030257.g.vcf \
  -O results/variants/SRR030257_filtered.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "basic_filter"

Explanation





VariantFiltration: Applies filters based on quality metrics:





QD < 2.0: Low quality per depth.



FS > 60.0: High strand bias.



MQ < 40.0: Low mapping quality.



Output: SRR030257_filtered.vcf.

Step 8: Variant Annotation

Annotate variants with functional information using Ensembl VEP.

Code

vep -i results/variants/SRR030257_filtered.vcf \
    -o results/variants/SRR030257_annotated.vcf \
    --cache --dir_cache /path/to/vep/cache \
    --assembly ASM584v2 \
    --everything

Explanation





VEP: Annotates variants with gene, transcript, and functional impact (e.g., missense, stop-gained).



Requires VEP cache for E. coli (ASM584v2). Install via ensembl-vep or use Ensembl’s online server for small datasets.



Output: SRR030257_annotated.vcf with biological context.

Step 9: Quality Control and Visualization

Aggregate QC metrics and visualize results.

Code

# Aggregate QC with MultiQC
multiqc results/ -o results/multiqc

# Generate coverage plot with mosdepth
mosdepth -n results/coverage/SRR030257 results/aligned/SRR030257_recal.bam

Explanation





MultiQC: Combines FastQC, Trimmomatic, and other logs into a single report (results/multiqc/multiqc_report.html).



mosdepth: Calculates read depth across the genome for coverage analysis.



For visualization, use IGV (manually load BAM/VCF) or R for plots (e.g., coverage histograms).

Step 10: Automation and Scalability

For production, automate with Nextflow (e.g., nf-core/sarek).

Code

# Example Nextflow command
nextflow run nf-core/sarek --input data/raw/*.fastq.gz --genome ASM584v2 -profile docker

Explanation





nf-core/sarek: A production-ready WGS pipeline with built-in QC, alignment, and variant calling.



Run on cloud platforms (e.g., AWS HealthOmics) for large datasets.



Configure with nextflow.config for custom references.

Conclusion

This pipeline processes WGS data from raw reads to annotated variants. Key outputs:





QC reports: results/fastqc/, results/multiqc/.



Aligned BAM: results/aligned/SRR030257_recal.bam.



Variants: results/variants/SRR030257_annotated.vcf.

For human WGS, use hg38 and public variant databases (e.g., dbSNP, 1000 Genomes). Adapt paths and parameters as needed. Push to GitHub (bioinformatics-pipelines) with a README.md detailing setup and rendering instructions.

References





Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data.



Bolger, A. M., et al. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.



Van der Auwera, G. A., et al. (2013). From FastQ data to high-confidence variant calls: The Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics, 43(1), 11.10.1-11.10.33.



McLaren, W., et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122.



Ewels, P., et al. (2020). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 36(3), 896-898.

(Note: Save references in references.bib for Quarto rendering.)
