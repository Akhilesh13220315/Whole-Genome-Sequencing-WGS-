# Whole Genome Sequencing Pipeline 

We'll use a small bacterial WGS dataset from NCBI SRA (SRR030257 - E. coli, ~200MB, suitable for Colab). This is open data from NCBI. For human WGS, datasets are larger, so bacterial is used for demo. Subsample if needed.

### Markdown Cell 1: Introduction
```
# End-to-End Whole Genome Sequencing Pipeline

This notebook implements an end-to-end WGS pipeline using open data from NCBI SRA.
- Dataset: SRR030257 (E. coli WGS, paired-end)
- Tools: FastQC, Trimmomatic, BWA, SAMtools, Picard, GATK, VEP (via pip/apt)
- Reference: E. coli genome (download from NCBI)

Run cells sequentially. Installations may take a few minutes.
```

### Code Cell 1: Install Dependencies
```
!apt-get update -qq
!apt-get install -y -qq fastqc trimmomatic bwa samtools default-jre
!pip install gatk4 # Note: GATK via pip is limited; for full, download jar
!wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
!unzip gatk-4.5.0.0.zip
!mv gatk-4.5.0.0/gatk /usr/local/bin/
!pip install ensembl-vep # For VEP, but needs configuration
!apt-get install -y tabix
```

### Code Cell 2: Download Open Data and Reference
```
# Download SRA data using fasterq-dump (install sra-toolkit)
!apt-get install -y sra-toolkit
!fasterq-dump SRR030257 --split-files -O data/

# Download E. coli reference (open from NCBI)
!mkdir ref
!wget -O ref/ecoli.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
!gunzip ref/ecoli.fa.gz
!bwa index ref/ecoli.fa
!samtools faidx ref/ecoli.fa
```

### Code Cell 3: Quality Control
```
!mkdir results
!fastqc -o results/fastqc data/SRR030257_1.fastq data/SRR030257_2.fastq
```

### Code Cell 4: Trimming
```
!trimmomatic PE -phred33 data/SRR030257_1.fastq data/SRR030257_2.fastq \
  results/trimmed/SRR030257_1_trimmed.fastq results/trimmed/SRR030257_1_unpaired.fastq \
  results/trimmed/SRR030257_2_trimmed.fastq results/trimmed/SRR030257_2_unpaired.fastq \
  ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Code Cell 5: Alignment
```
!bwa mem -R "@RG\tID:SRR030257\tSM:SRR030257\tLB:lib1\tPL:illumina" \
  ref/ecoli.fa results/trimmed/SRR030257_1_trimmed.fastq results/trimmed/SRR030257_2_trimmed.fastq > results/aligned/SRR030257.sam
!samtools view -bS results/aligned/SRR030257.sam | samtools sort -o results/aligned/SRR030257_sorted.bam
!samtools index results/aligned/SRR030257_sorted.bam
```

### Code Cell 6: Post-Alignment (Mark Duplicates, BQSR)
```
# Picard for duplicates (install Picard)
!wget https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar
!java -jar picard.jar MarkDuplicates I=results/aligned/SRR030257_sorted.bam O=results/aligned/SRR030257_dedup.bam M=results/aligned/metrics.txt

# GATK BQSR (needs known sites; for demo, skip or use dummy)
!gatk BaseRecalibrator -I results/aligned/SRR030257_dedup.bam -R ref/ecoli.fa --known-sites known.vcf -O recal.table # Assume known.vcf exists
!gatk ApplyBQSR -I results/aligned/SRR030257_dedup.bam -R ref/ecoli.fa --bqsr-recal-file recal.table -O results/aligned/SRR030257_recal.bam
```

### Code Cell 7: Variant Calling
```
!gatk HaplotypeCaller -R ref/ecoli.fa -I results/aligned/SRR030257_recal.bam -O results/variants/SRR030257.vcf
```

### Code Cell 8: Annotation (VEP setup needed)
```
!vep -i results/variants/SRR030257.vcf -o results/variants/annotated.vcf --cache --assembly ASM584v2
```

### Markdown Cell 2: Conclusion
```
Pipeline complete. Download results from /results folder.
```

Upload to GitHub repo as `wgs_colab.ipynb`.

# RNA-Seq Pipeline in Google Colab

Use open data from GEO GSE18508 (small subsampled RNA-Seq for D. melanogaster DE). Files from Zenodo (open).

### Markdown Cell 1: Introduction
```
# End-to-End RNA-Seq Pipeline

Using subsampled open data from GEO GSE18508 for differential expression.
- Tools: FastQC, Trimmomatic, STAR or Salmon, featureCounts, DESeq2
- Reference: D. melanogaster genome
```

### Code Cell 1: Install Dependencies
```
!apt-get update -qq
!apt-get install -y fastqc trimmomatic hisat2 samtools subread r-base
!pip install multiqc HTSeq salmon
!R -e "install.packages(c('DESeq2', 'pheatmap'), repos='http://cran.rstudio.com/')"
```

### Code Cell 2: Download Data
```
!mkdir data
!wget -P data/ https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger
!wget -P data/ https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger
!wget -P data/ https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger
!wget -P data/ https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger

# Download reference (D. melanogaster)
!mkdir ref
!wget -O ref/dmel.fa.gz ftp://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
!gunzip ref/dmel.fa.gz
!wget -O ref/dmel.gtf https://zenodo.org/record/6457007/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
!gunzip ref/dmel.gtf
!hisat2-build ref/dmel.fa ref/dmel_index # For HISAT2 alternative if STAR not installed
```

### Code Cell 3: QC and Trimming
```
!fastqc -o results/fastqc data/*.fastqsanger
!for sample in GSM461177 GSM461180; do
  trimmomatic PE data/${sample}_1_subsampled.fastqsanger data/${sample}_2_subsampled.fastqsanger \
    results/trimmed/${sample}_1_trim.fastq results/trimmed/${sample}_1_unp.fastq \
    results/trimmed/${sample}_2_trim.fastq results/trimmed/${sample}_2_unp.fastq \
    ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 MINLEN:25
done
```

### Code Cell 4: Quantification with Salmon
```
!salmon index -t ref/dmel.fa -i ref/salmon_index
!for sample in GSM461177 GSM461180; do
  salmon quant -i ref/salmon_index -l A -1 results/trimmed/${sample}_1_trim.fastq -2 results/trimmed/${sample}_2_trim.fastq -o results/quant/${sample}
done
```

### Code Cell 5: Differential Expression with DESeq2
```
%%R
library(DESeq2)
# Assume count files from salmon or featureCounts
# For simplicity, manual counts or use tximport
dir <- "results/quant"
samples <- c("GSM461177", "GSM461180")
condition <- c("untreated", "treated")
colData <- data.frame(sample=samples, condition=condition)
dds <- DESeqDataSetFromTximport(tximport(dir, type="salmon"), colData, ~condition) # Assume tximport installed
dds <- DESeq(dds)
res <- results(dds)
head(res)
```

Upload as `rna_seq_colab.ipynb` to GitHub.

# ChIP-Seq Pipeline in Google Colab

Use open data from ENCODE (small subsets if possible; here use accessions and subsample with seqtk).

From browse, accessions like ENCFF000PED for HeLa ChIP.

### Markdown Cell 1: Introduction
```
# End-to-End ChIP-Seq Pipeline

Using open ChIP-Seq data from ENCODE for CTCF.
- Dataset: ENCFF000PED (HeLa ChIP), ENCFF000PET (input)
- Tools: FastQC, Bowtie2, MACS2
```

### Code Cell 1: Install Dependencies
```
!apt-get update -qq
!apt-get install -y fastqc bowtie2 samtools macs2 seqtk
!pip install multiqc
```

### Code Cell 2: Download Data
```
!mkdir data
!wget -O data/chip.fastq.gz https://www.encodeproject.org/files/ENCFF000PED/@@download/ENCFF000PED.fastq.gz
!wget -O data/input.fastq.gz https://www.encodeproject.org/files/ENCFF000PET/@@download/ENCFF000PET.fastq.gz

# Subsample to 1M reads for small analysis
!seqtk sample -s100 data/chip.fastq.gz 1000000 > data/chip_small.fastq
!seqtk sample -s100 data/input.fastq.gz 1000000 > data/input_small.fastq
```

### Code Cell 3: Alignment (use hg38 reference)
```
!mkdir ref
!wget -O ref/hg38.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
!gunzip ref/hg38.fa.gz
!bowtie2-build ref/hg38.fa ref/hg38_index
!bowtie2 -x ref/hg38_index -U data/chip_small.fastq -S results/chip.sam
!samtools view -bS results/chip.sam | samtools sort -o results/chip.bam
!bowtie2 -x ref/hg38_index -U data/input_small.fastq -S results/input.sam
!samtools view -bS results/input.sam | samtools sort -o results/input.bam
!samtools index results/chip.bam
!samtools index results/input.bam
```

### Code Cell 4: Peak Calling
```
!macs2 callpeak -t results/chip.bam -c results/input.bam -f BAM -g hs --outdir results/peaks -n sample
```

Upload as `chip_seq_colab.ipynb`.

# ATAC-Seq Pipeline in Google Colab

Use open data from ENCODE (ENCFF398QLV for NK cells).

### Markdown Cell 1: Introduction
```
# End-to-End ATAC-Seq Pipeline

Using open ATAC-Seq data from ENCODE.
- Dataset: ENCFF398QLV (NK cell untreated)
- Tools: FastQC, Bowtie2, MACS2, deepTools
```

### Code Cell 1: Install Dependencies
```
!apt-get update -qq
!apt-get install -y fastqc bowtie2 samtools macs2
!pip install deeptools
```

### Code Cell 2: Download Data
```
!mkdir data
!wget -O data/atac1.fastq.gz https://www.encodeproject.org/files/ENCFF398QLV/@@download/ENCFF398QLV.fastq.gz
!wget -O data/atac2.fastq.gz https://www.encodeproject.org/files/ENCFF363HBZ/@@download/ENCFF363HBZ.fastq.gz

# Subsample
!seqtk sample -s100 data/atac1.fastq.gz 1000000 > data/atac1_small.fastq
!seqtk sample -s100 data/atac2.fastq.gz 1000000 > data/atac2_small.fastq
```

### Code Cell 3: Alignment
```
# Use same hg38 ref as above
!bowtie2 -x ref/hg38_index -U data/atac1_small.fastq -S results/atac1.sam
!samtools view -bS results/atac1.sam | samtools sort -o results/atac1.bam
!bowtie2 -x ref/hg38_index -U data/atac2_small.fastq -S results/atac2.sam
!samtools view -bS results/atac2.sam | samtools sort -o results/atac2.bam
```

### Code Cell 4: Filtering and Peak Calling
```
# ATAC shift
!alignmentSieve --ATACshift --bam results/atac1.bam -o results/atac1_shifted.bam
!macs2 callpeak -t results/atac1_shifted.bam -f BAM --nomodel --shift -75 --extsize 150 -g hs -n atac1 --outdir results/peaks
```

Upload as `atac_seq_colab.ipynb`.

To add to GitHub repo, save these as .ipynb files and push to `bioinformatics-pipelines`. Use open data links as shown. For large files, subsampling keeps it feasible in Colab. If errors, adjust paths or install missing tools.