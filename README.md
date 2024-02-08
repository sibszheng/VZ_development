# ChIP-seq pipeline (using TFEB_CHIP_rep1 as an example)
## quality control
fastp -i TFEB_ChIP_rep1_R1.fastq.gz \
      -I TFEB_ChIP_rep1_R2.fastq.gz \
      -o TFEB_ChIP_rep1_R1.fastp.fastq.gz \
      -O TFEB_ChIP_rep1_R2.fastp.fastq.gz \
      -w 8 -l 25
## genome mapping
hisat2 -p 6 \
          -mm10_index \
          --no-temp-splicesite \
          --no-spliced-alignment \
          --summary-file TFEB_ChIP_rep1_mapping_stats.txt \
          -1 TFEB_ChIP_rep1_R1.fastp.fastq.gz \
          -2 TFEB_ChIP_rep1_R2.fastp.fastq.gz | \
          samtools view -@ 6 -ShuF 4 -q 30 - | \
          samtools sort TFEB_ChIP_rep1_tmp -o TFEB_ChIP_rep1_sorted.bam
samtools index TFEB_ChIP_rep1_sorted.bam
## duplicate removal
java -Xmx5g \
       -XX:ParallelGCThreads=8 \
       -jar picard.jar MarkDuplicates \
       I= Spm_H3K4me3_rep1.sorted.bam \
       O= Spm_H3K4me3_rep1.sorted.picardMD.bam \
       M= Spm_H3K4me3_rep1.sorted.picardMD.txt \
       REMOVE_DUPLICATES=true 
