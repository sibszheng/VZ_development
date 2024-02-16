#There are three replicates for each condition(TFEB_NC_rep1,2,3 and TFEB_KD_rep1,2,3). Here TFEB_KD_rep1 is used as an example.
#1 Quality control
fastp -i TFEB_KD_rep1_R1.fastq.gz \
-I TFEB_KD_rep1_R2.fastq.gz \
-o TFEB_KD_rep1_R1.fastp.fastq.gz \
-O TFEB_KD_rep1_R2.fastp.fastq.gz \
-w 8 -l 25

#2 Genome mapping
hisat2 -p 6 \
-mm10_index \
--summary-file TFEB_KD_rep1_mapping_stats.txt \
-1 TFEB_KD_rep1_R1.fastp.fastq.gz \
-2 TFEB_KD_rep1_R2.fastp.fastq.gz | \
samtools view -@ 6 -ShuF 4 -q 30 - | \
samtools sort TFEB_KD_rep1_tmp -o TFEB_KD_rep1_sorted.bam
samtools index TFEB_KD_rep1_sorted.bam

#3 Duplicate removal (picard)
java -jar -Xmx8g picard.jar MarkDuplicates \
INPUT= TFEB_KD_rep1_sorted.bam \
OUPUT= TFEB_KD_rep1_sorted_picardMD.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE= TFEB_KD_rep1_sorted_picardMD.met
samtools index TFEB_KD_rep1_sorted_picardMD.bam

#4 Differential gene expression analysis
#Gene count matrix was obtained using htseq then differential gene expression analysis was performed using DESeq2 in R.
htseq-count -s no -f bam -i gene_name \
TFEB_KD_rep1_sorted_picardMD.bam \
gencode.vM10.annotation.gtf \ 
> TFEB_KD_rep1_count
paste TFEB_NC_reps_count TFEB_KD_reps_count | \
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
res_TFEB_KD <- results(dds,contrast = c("treatment","KD","NC"))
TFEB_KD_up <- as.data.frame(subset(res_TFEB_KD,padj<0.001 & log2FoldChange>0.5))
TFEB_KD_down <- as.data.frame(subset(res_TFEB_KD,padj<0.001 & log2FoldChange<-0.5))

#5 GO enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db)
TFEB_KD_up_BP <- enrichGO(gene=rownames(TFEB_KD_up), OrgDb=org.Mm.eg.db,ont='BP',keyType = 'SYMBOL')

#6 Calculation for FPKM
#Gene count and length was obtained using featureCount then subjected to calculation for FPKM in R.
featureCounts -F gtf -t exon -g gene_name -p -T 4 \
TFEB_KD_rep1_sorted_picardMD.bam \
-a gencode.vM10.annotation.gtf \
-o TFEB_KD_rep1_count_length 

KD_rep1_matrix <- read.table (file = 'TFEB_KD_rep1_count_length', header = T)
colnames(KD_rep1_matrix)[,7] <- "count¡±
scaling_factor <- colSums(KD_rep1_matrix$count)/1000000
KD_rep1_matrix$FPKM <- (1000 * KD_rep1_matrix$count) / (scaling_factor * KD_rep1_matrix $Length.x)

