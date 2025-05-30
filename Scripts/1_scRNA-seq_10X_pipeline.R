#There are four biological replicates:P0_CD133_RGC_Droplet, P5_CD133_RGC_Droplet, P0_CD133_RGC_Droplet_rep and P5_CD133_RGC_Droplet_rep.Here P0_CD133_RGC_Droplet_rep file was used as an example.

#1 Building the genome index of mm10
STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir STAR_index/mm10/ \
--genomeFastaFiles mm10.fa \
--sjdbGTFfile gencode.vM10.annotation.gtf
#The gtf file can be found in the Mapping folder.

#2 Genome mapping
STAR --runThreadN 4 \
--genomeDir STAR_index/mm10/ \
--readFilesCommand zcat \
--outFileNamePrefix STAR_outs/ \
--readFilesIn P0_CD133_RGC_Droplet_rep_Read2 P0_CD133_RGC_Droplet_rep_Read1 \
--soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
--soloCBwhitelist whitelist_10X.txt \
--soloCellFilter EmptyDrops_CR \
--soloStrand Forward \
--outSAMattributes CB UB \
--outSAMtype BAM SortedByCoordinate \
--soloBarcodeReadLength 0
#The whitellist file can be downloaded from https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz.
#In the output ¡°Solo.out¡± directory, there are ¡°Gene/raw¡± and ¡°Gene/filtered¡± folders, which contain expression values of each gene in each cell. The files in the ¡°Gene/raw¡± directory contain the gene expression of every barcode in the whitelist. The files in the ¡°Gene/filtered¡± contain similar information with barcodes that have too few reads removed. The files in the ¡°Gene/filtered¡± were imported into R for downstream analysis.

#3 Downstream analysis in R
#3-1 Loading the required R packages
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(patchwork)
library(SingleR)
library(monocle)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

#3-2 Doublet removal (Seurat and Doubletfinder)
#3-2-1 Creating Seurat object
P0_counts_rep2 <- Read10X(data.dir = '/Solo.out/Gene/filtered')
sce_P0_rep2 <- CreateSeuratObject(counts = P0_counts_rep2)
#3-2-2 Normalization
sce_P0_rep2 <- NormalizeData(sce_P0_rep2) 
sce_P0_rep2 <- FindVariableFeatures(sce_P0_rep2, selection.method = "vst",nfeatures = 2000)
sce_P0_rep2 <- ScaleData(sce_P0_rep2,vars.to.regress ="nCount_RNA")
#3-2-3 Linear dimension reduction
sce_P0_rep2 <- RunPCA(sce_P0_rep2,ndims.print = 1:5, nfeatures.print = 5)
#3-2-4 Clustering
sce_P0_rep2 <- FindNeighbors(sce_P0_rep2, dims = 1:20)
sce_P0_rep2 <- FindClusters(sce_P0_rep2, resolution = 0.8)
sce_P0_rep2 <- RunTSNE(object = sce_P0_rep2, dims.use = 1:20)
DimPlot(sce_P0_rep2,reduction="tsne",label=T)
#3-2-5 Doublet identification
sweep.res.list <- paramSweep_v3(sce_P0_rep2, PCs = 1:20)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = 0.076
homotypic.prop <- modelHomotypic(sce_P0_rep2$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(sce_P0_rep2)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
sce_P0_rep2 <- doubletFinder_v3(sce_P0_rep2, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F)
DimPlot(sce_P0_rep2, reduction = "tsne", group.by = "DF.classifications_0.25_0.005_584")
#3-2-6 Doublet removal
sce_P0_rep2 <- CreateSeuratObject(counts=sce_P0_rep2[,sce_P0_rep2$DF.classifications_0.25_0.005_584=="Singlet"]@assays[["RNA"]]@counts, project = "10X_P0_rep2")

#3-3 
sce_P0_rep2 <- NormalizeData(sce_P0_rep2)
sce_P0_rep2 <- FindVariableFeatures(sce_P0_rep2, selection.method = "vst",nfeatures = 2000)
sce_P0_rep2 <- ScaleData(sce_P0_rep2,vars.to.regress ="nCount_RNA") 
sce_P0_rep2 <- RunPCA(sce_P0_rep2,ndims.print = 1:5, nfeatures.print = 5)
sce_P0_rep2 <- FindNeighbors(sce_P0_rep2, dims = 1:20)
sce_P0_rep2 <- FindClusters(sce_P0_rep2, resolution = 0.8)
sce_P0_rep2 <- RunTSNE(object = sce_P0_rep2, dims.use = 1:20)
DimPlot(sce_P0_rep2,reduction="tsne",label=T)

#3-4 Celltype annotation (SingleR)
counts <- GetAssayData(sce_P0_rep2)
mg <- celldex::MouseRNAseqData()
common_genes <- intersect(rownames(counts), rownames(mg)) 
counts <- counts[common_genes,]
mg <- mg[common_genes,]
celltype <- SingleR(test = counts, ref = mg, labels = mg$label.fine) 
sce_P0_rep2$celltype <- celltype$labels

#3-5 Calculation for mitochondrial gene percentage
MITO <- read.csv("Mouse.Mitocarta3.0.csv",header = FALSE)
counts <- as.matrix(sce_P0_rep2@assays[["RNA"]]@counts)
mitopercent <- apply(counts,2,function(x){sum(x>0&rownames(counts)[x]%in%MITO[,1])}) /sce_P0_rep2$nFeature_RNA*100
sce_P0_rep2$mitopercent <- mitopercent

#3-6 Filtering endothelial cells, microglia, pericytes/fibroblasts and dying cells (Monocle and Seurat)
#3-6-1 Creating Monocle object
counts <- sce_P0_rep2@assays[["RNA"]]@counts 
gene_ann <- data.frame(gene_short_name=row.names(counts),row.names = row.names(counts))
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd <- new("AnnotatedDataFrame",data= sce_P0_rep2@meta.data)
cds_P0_rep2<- newCellDataSet(counts,phenoData=pd,featureData=fd,
                        expressionFamily=negbinomial.size())
#3-6-2 Normalization
cds_P0_rep2<- estimateSizeFactors(cds_P0_rep2)
cds_P0_rep2<- estimateDispersions(cds_P0_rep2)
#3-6-3 Clustering
disp_table <- dispersionTable(cds_P0_rep2)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_P0_rep2<- setOrderingFilter(cds_P0_rep2, unsup_clustering_genes$gene_id)
cds_P0_rep2<- reduceDimension(cds_P0_rep2, max_components = 2, num_dim = 6,
                           reduction_method = 'tSNE', verbose = T)
cds_P0_rep2<- clusterCells(cds_P0_rep2, num_clusters = 20)
plot_cell_clusters(cds_P0_rep2)
#3-6-4 Filtering
table(cds_P0_rep2X$Cluster, cds_P0_rep2$celltype) #Cluster5,6 is constituted of endothelial cells, microglia and pericytes/fibroblasts, which have been annotated by SingleR.
sce_P0_rep2$monocle_clusters <- cds_P0_rep2$Cluster
VlnPlot(sce_P0_rep2,features=c("nFeature_RNA","nCount_RNA","mitopercent"),group.by="monocle_clusters")
sce_P0_rep2<- sce_P0_rep2[,!(sce_P0_rep2$monocle_clusters=="5"|
                       sce_P0_rep2$monocle_clusters=="6")]
#Cells of Cluster5 and Cluster6 are excluded from subsequent analysis.

#3-7 Integrating rep1 and rep2 data (Seurat)
sce_10X_all <- merge(sce_P0_rep2, y = c(sce_P0_rep1,sce_P5_rep1,sce_P5_rep2),project = "sce_10X_all") # All data had been filtered.
sce_10X_integrated <- NormalizeData(sce_10X_all)
sce_10X_integrated <- FindVariableFeatures(sce_10X_integrated)
sce_10X_integrated <- ScaleData(sce_10X_integrated)
sce_10X_integrated <- RunPCA(sce_10X_integrated)
sce_10X_integrated <- IntegrateLayers(object = sce_10X_integrated, method=HarmonyIntegration,orig.reduction = "pca",new.reduction ="harmony",verbose = FALSE)
sce_10X_integrated <- RunTSNE(sce_10X_integrated, reduction = "harmony", dims =1:30, reduction.name = "tsne.harmony")

#3-8 Clustering (Seurat)
sce_10X_integrated <- FindNeighbors(sce_10X_integrated, reduction = "harmony", dims = 1:30)
sce_10X_integrated <- FindClusters(sce_10X_integrated, resolution = 0.4, cluster.name = "harmony_clusters")

#3-9 Trajectory analysis (Monocle)
counts <- sce_10X_all@assays[["RNA"]]@counts 
gene_ann <- data.frame(gene_short_name=row.names(counts),row.names = row.names(counts))
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd <- new("AnnotatedDataFrame",data= sce_10X_integrated@meta.data)
cds_10X_integrated <- newCellDataSet(counts,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
cds_10X_integrated <- estimateSizeFactors(cds_10X_integrated)
cds_10X_integrated <- estimateDispersions(cds_10X_integrated)
diff_test_res <- differentialGeneTest(cds_10X_integrated,fullModelFormulaStr ="~harmony_clusters",cores = 10)
ordering_genes <- diff_test_res$gene_short_name[order(diff_test_res$qval)][1:2000]
cds_10X_integrated <- setOrderingFilter(cds_10X_integrated, ordering_genes)
cds_10X_integrated <- reduceDimension(cds_10X_integrated, max_components = 2, method='DDRTree') 
cds_10X_integrated <- orderCells(cds_10X_integrated)
cds_10X_integrated <- orderCells(cds_10X_integrated,root_state = 2)#State1,2,3 represent nEPC state, bGPC state and nNSC-NB state, respectively.
plot_cell_trajectory(cds_10X_integrated)

#3-10 Differential expression analysis on branches (Monocle)
Branch_diff <- BEAM(cds_10X_integrated)
Branch_diff<- Branch_diff[order(Branch_diff$qval),]

#3-11 Differential expression analysis on clusters (Seurat)
sce_10X_all <- NormalizeData(sce_10X_all) 
sce_10X_all <- FindVariableFeatures(sce_10X_all, nfeatures = 2000)
sce_10X_all <- ScaleData(sce_10X_all,vars.to.regress ="nCount_RNA") 
sce_10X_all@harmonny_clusters<-sce_10X_integrated$harmonny_clusters
#3-11-1 Marker gene identification for NSC clusters
sce_10X_NSC<-sce_10X_all[,sce_10X_all$harmonny_clusters ==2|
                          sce_10X_all$harmonny_clusters ==4|
                          sce_10X_all$harmonny_clusters ==6|
                          sce_10X_all$harmonny_clusters ==12]
NSC_cluster_markers<-FindAllMarkers(sce_10X_NSC,only.pos = TRUE)
#3-11-2 Marker gene identification for EPC clusters
sce_10X_EPC<-sce_10X[,sce_10X_all$harmonny_clusters ==8|
                      sce_10X_all$harmonny_clusters ==9]
EPC_cluster_markers<-FindAllMarkers(sce_10X_EPC,only.pos = TRUE)

#3-12 GO enrichment analysis (clusterProfiler)
nNSC_early_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1 
                                               & NSC_cluster_markers$cluster==4]
nNSC_late_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1
                                              & NSC_cluster_markers$cluster==6]
NB1_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1
                                        & NSC_cluster_markers$cluster==12]
NB2_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1
                                        & NSC_cluster_markers$cluster==2]
nEPC_early_markers <- EPC_cluster_markers$gene[EPC_cluster_markers$avg_log2FC>1
                                              & EPC_cluster_markers$cluster==8]
nEPC_late_markers <- EPC_cluster_markers$gene[EPC_cluster_markers$avg_log2FC>1
                                             & EPC_cluster_markers$cluster==9]
nNSC_early_markers_BP <- enrichGO(gene=nNSC_early_markers,
                                 OrgDb=org.Mm.eg.db,ont='BP',
                                 keyType = 'SYMBOL')
nNSC_late_markers_BP <- enrichGO(gene=nNSC_late_markers,
                              OrgDb=org.Mm.eg.db,ont='BP',
                              keyType = 'SYMBOL')
NB1_markers_BP <- enrichGO(gene=NB1_markers,
                         OrgDb=org.Mm.eg.db,ont='BP',
                         keyType= 'SYMBOL')
NB2_markers_BP <- enrichGO(gene=NB2_markers,
                         OrgDb=org.Mm.eg.db,ont='BP',
                         keyType= 'SYMBOL')
nEPC_early_markers_BP <- enrichGO(gene=nEPC_early_markers,
                               OrgDb=org.Mm.eg.db,ont='BP',
                               keyType = 'SYMBOL')
nEPC_late_markers_BP <- enrichGO(gene=nEPC_late_markers,
                              OrgDb=org.Mm.eg.db,ont='BP',
                              keyType = 'SYMBOL')
