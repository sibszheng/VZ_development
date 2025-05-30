#There are three replicates for each condition. Here siTfeb_rep1 is used as an example.
#1 Quality control
fastp -i siTfeb_rep1_R1.fastq.gz \
-I siTfeb_rep1_R2.fastq.gz \
-o siTfeb_rep1_R1.fastp.fastq.gz \
-O siTfeb_rep1_R2.fastp.fastq.gz \
-w 8 -l 25

#2 Genome mapping
hisat2 -p 6 \
-mm10_index \
--summary-file siTfeb_rep1_mapping_stats.txt \
-1 siTfeb_rep1_R1.fastp.fastq.gz \
-2 siTfeb_rep1_R2.fastp.fastq.gz | \
samtools view -@ 6 -ShuF 4 -q 30 - | \
samtools sort siTfeb_rep1_tmp -o siTfeb_rep1_sorted.bam
samtools index siTfeb_rep1_sorted.bam

#3 Duplicate removal (picard)
java -jar -Xmx8g picard.jar MarkDuplicates \
INPUT= siTfeb_rep1_sorted.bam \
OUPUT= siTfeb_rep1_sorted_picardMD.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE= siTfeb_rep1_sorted_picardMD.met
samtools index siTfeb_rep1_sorted_picardMD.bam

#4 Differential gene expression analysis
#Gene count matrix was obtained using htseq then differential gene expression analysis was performed using DESeq2 in R.
htseq-count -s no -f bam -i gene_name \
siTfeb_rep1_sorted_picardMD.bam \
gencode.vM10.annotation.gtf \ 
> siTfeb_rep1_count
paste siNC_reps_count siTfeb_reps_count | \
cut -f 1-2,4,6,8,10,12 | \
head -48321 < TFEB_RNAseq_count

library(DESeq2)
cts <- read.csv("TFEB_RNAseq_count",sep="\t",header = FALSE,
                col.names = c("gene_name", "NC1","NC2","NC3","KD1","KD2","KD3"))
rownames(cts) <- cts[,1]
cts <- cts[,-1]
condition <- c("NC","NC","NC","KD","KD","KD")
condition <- as.data.frame(condition)
rownames(condition) <- colnames(cts)
colnames(condition) <- "treatment"
condition$treatment <- factor(condition$treatment)
cts <- cts[rowSums(cts)>1,]
dds <- DESeqDataSetFromMatrix(countData=cts,colData=condition,design = ~treatment)
dds<-DESeq(dds)
ntd <- normTransform(dds) 
res_siTfeb <- results(dds,contrast = c("treatment","KD","NC"))
siTfeb_up <- as.data.frame(subset(res_siTfeb,padj<0.001 & log2FoldChange>0.5))
siTfeb_down <- as.data.frame(subset(res_siTfeb,padj<0.001 & log2FoldChange<-0.5))

#5 GO enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db)
siTfeb_up_BP <- enrichGO(gene=rownames(siTfeb_up), OrgDb=org.Mm.eg.db,ont='BP',keyType = 'SYMBOL')

#6 Calculation for FPKM
#Gene count and length was obtained using featureCount then subjected to calculation for FPKM in R.
featureCounts -F gtf -t exon -g gene_name -p -T 4 \
siTfeb_rep1_sorted_picardMD.bam \
-a gencode.vM10.annotation.gtf \
-o siTfeb_rep1_count_length 

KD_rep1_matrix <- read.table (file = 'siTfeb_rep1_count_length', header = T)
colnames(KD_rep1_matrix)[,7] <- "count¡±
scaling_factor <- colSums(KD_rep1_matrix$count)/1000000
KD_rep1_matrix$FPKM <- (1000 * KD_rep1_matrix$count) / (scaling_factor * KD_rep1_matrix $Length.x)

