# Sample data set from encode

## Download reference:

Download the reference to map the fastq files:

wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

conda create --name bwa_samtools
conda activate bwa_samtools
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bedtools

gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

bwa index Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

* Adapted from [here](https://www.biostars.org/p/369483/) 

* Some info got from the following Biostar question:
[Example datasets for RNA-Seq Differential Expression Analysis](https://www.biostars.org/p/368224/)

```bash
conda activate bwa_samtools
bwa mem -M -Y Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa SRR5093452_GSM2421617_polyA_mRNA_RNA-seq_from_A549_ENCLB575NVN_Homo_sapiens_RNA-Seq_1.fastq.gz SRR5093452_GSM2421617_polyA_mRNA_RNA-seq_from_A549_ENCLB575NVN_Homo_sapiens_RNA-Seq_2.fastq.gz |  
samtools view -S -b > SRR5093452.bam

# [main] Real time: 6491.282 sec; CPU: 6573.601 sec

samtools sort -n -O bam SRR5093452.bam -o SRR5093452_sorted.bam

samtools index SRR5093452_sorted.bam

samtools view -b SRR5093452_sorted.bam chr22 | samtools fastq -1 SRR5093452_chr22_1.fastq.gz -2 SRR5093452_chr22_2.fastq.gz
```

ftp://ftp.broadinstitute.org/pub/seq/references/

bwa index Homo_sapiens_assembly19.fasta

bwa mem Homo_sapiens_assembly19.fasta SRR5093452_GSM2421617_polyA_mRNA_RNA-seq_from_A549_ENCLB575NVN_Homo_sapiens_RNA-Seq_1.fastq.gz SRR5093452_GSM2421617_polyA_mRNA_RNA-seq_from_A549_ENCLB575NVN_Homo_sapiens_RNA-Seq_2.fastq.gz |  
samtools sort -o SRR5093452.bam

# already sorted in previous step 
samtools index SRR5093452.bam

samtools view  -b SRR5093452.bam 22 > SRR5093452_chr22.bam 
samtools fastq -1 SRR5093452_chr22_1.fastq.gz -2 SRR5093452_chr22_2.fastq.gz SRR5093452_chr22.bam


## FASTA chromosome 22

samtools faidx Homo_sapiens_assembly19.fasta
samtools faidx Homo_sapiens_assembly19.fasta 22  > Homo_sapiens_assembly19_chr22.fasta