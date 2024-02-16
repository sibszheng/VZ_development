#There are two fastq files:VZ_P0_10X and VZ_P5_10X.Here VZ_P0_10X file was used as an example.

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
--readFilesIn VZ_P0_10X_R2.fastq.gz VZ_P0_10X_R1.fastq.gz \
--soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
--soloCBwhitelist whitelist_10X.txt \
--soloCellFilter EmptyDrops_CR \
--soloStrand Forward \
--outSAMattributes CB UB \
--outSAMtype BAM SortedByCoordinate \
--soloBarcodeReadLength 0
#The whitellist file can be downloaded from https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz.
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
library(destiny)

#3-2 Doublet removal (Seurat and Doubletfinder)
#3-2-1 Creating Seurat object
P0_counts <- Read10X(data.dir = '/Solo.out/Gene/filtered')
sce_P0 <- CreateSeuratObject(counts = P0_counts, project = "sce_P0")
#3-2-2 Normalization
sce_P0 <- NormalizeData(sce_P0) 
sce_P0 <- ScaleData(sce_P0,vars.to.regress ="nCount_RNA")
#3-2-3 Linear dimension reduction
sce_P0 <- FindVariableFeatures(sce_P0, selection.method = "vst",nfeatures = 2000)
sce_P0 <- RunPCA(sce_P0,ndims.print = 1:5, nfeatures.print = 5)
#3-2-4 Clustering
sce_P0 <- FindNeighbors(sce_P0, dims = 1:20)
sce_P0 <- FindClusters(sce_P0, resolution = 0.8)
sce_P0 <- RunTSNE(object = sce_P0, dims.use = 1:20)
DimPlot(sce_P0,reduction="tsne",label=T)
#3-2-5 Doublet identification
sweep.res.list <- paramSweep_v3(sce_P0, PCs = 1:20)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = 0.065
homotypic.prop <- modelHomotypic(sce_P0$monocle_clusters)
nExp_poi <- round(DoubletRate*ncol(sce_P0)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
sce_P0 <- doubletFinder_v3(sce_P0, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                           nExp = nExp_poi.adj, reuse.pANN = F)
DimPlot(sce_P0, reduction = "tsne", group.by = "DF.classifications_0.25_0.11_512")
#3-2-6 Doublet removal
sce.singlet_P0 <- CreateSeuratObject(counts=sce_P0[,sce_P0$DF.classifications_0.25_0.11_512=="Singlet"]@assays[["RNA"]]@counts, project = "P0_exp")

#3-3 Merging P0 and P5 data (Seurat)
sce_10X <- merge(sce.singlet_P0, y = sce.singlet_P5, 
                 add.cell.ids = c("P0", "P5"), project = "sce_10X")
sce_10X <- NormalizeData(sce_10X) 
sce_10X <- ScaleData(sce_10X,vars.to.regress ="nCount_RNA") 
sce_10X <- FindVariableFeatures(sce_10X, selection.method = "vst",nfeatures = 2000)
sce_10X <- RunPCA(sce_10X,ndims.print = 1:5, nfeatures.print = 5)
sce_10X <- FindNeighbors(sce_10X, dims = 1:20)
sce_10X <- FindClusters(sce_10X, resolution = 0.8)
sce_10X <- RunTSNE(object = sce_10X, dims.use = 1:20)
DimPlot(sce_10X,reduction="tsne",label=T)

#3-4 Celltype annotation (SingleR)
counts <- GetAssayData(sce_10X)
mg <- celldex::MouseRNAseqData()
common_genes <- intersect(rownames(counts), rownames(mg)) 
counts <- counts[common_genes,]
mg <- mg[common_genes,]
celltype <- SingleR(test = counts, ref = mg, labels = mg$label.fine) 
sce_10X$celltype <- celltype$labels

#3-5 Calculation for mitochondrial gene percentage
MITO <- read.csv("Mouse.Mitocarta3.0.csv",header = FALSE)
counts <- as.matrix(sce_10X@assays[["RNA"]]@counts)
mitopercent <- apply(counts,2,function(x){sum(x>0&rownames(bb)[x]%in%MITO[,1])}) /sce_10X$nFeature_RNA*100
sce_10X$mitopercent <- mitopercent

#3-6 Filtering endothelial cells, microglia, pericytes/fibroblasts and dying cells (Monocle and Seurat)
#3-6-1 Creating Monocle object
counts <- sce_10X@assays[["RNA"]]@counts 
gene_ann <- data.frame(gene_short_name=row.names(counts),row.names = row.names(counts))
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd <- new("AnnotatedDataFrame",data= sce_10X@meta.data)
cds_10X <- newCellDataSet(counts,phenoData=pd,featureData=fd,
                        expressionFamily=negbinomial.size())
#3-6-2 Normalization
cds_10X <- estimateSizeFactors(cds_10X)
cds_10X <- estimateDispersions(cds_10X)
#3-6-3 Clustering
disp_table <- dispersionTable(cds_10X)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_10X <- setOrderingFilter(cds_10X, unsup_clustering_genes$gene_id)
cds_10X <- reduceDimension(cds_10X, max_components = 2, num_dim = 6,
                           reduction_method = 'tSNE', verbose = T)
cds_10X <- clusterCells(cds_10X, num_clusters = 15)
plot_cell_clusters(cds_10X)
#3-6-4 Filtering
table(cds_10X$Cluster, cds_10X$celltype) #Cluster8 is constituted of endothelial cells, microglia and pericytes/fibroblasts, which have been annotated by SingleR.
sce_10X$monocle_clusters <- cds_10X$Cluster
VlnPlot(sce_10X,features=c("nFeature_RNA","nCount_RNA","mitopercent"),group.by="monocle_clusters")#Cells of Cluster3 show low gene number and UMI number but high mitochondrial gene percentage, indicating that they are dying cells.
cds_10X <- cds_10X [,!(cds_10X$Cluster=="8"| 
                       cds_10X$Cluster=="3")]
sce_10X<- sce_10X[,!(sce_10X$monocle_clusters=="8"|
                     sce_10X$monocle_clusters=="3")]
#Cells of Cluster8 and Cluster3 are excluded from subsequent analysis.

#3-7 Second-round filtering of dying cells (Monocle and Seurat)
cds_10X <- estimateSizeFactors(cds_10X)
cds_10X <- estimateDispersions(cds_10X)
disp_table <- dispersionTable(cds_10X)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_10X <- setOrderingFilter(cds_10X, unsup_clustering_genes$gene_id)
cds_10X <- reduceDimension(cds_10X, max_components = 2, num_dim = 6,
                           reduction_method = 'tSNE', verbose = T)
cds_10X <- clusterCells(cds_10X, num_clusters = 13)
plot_cell_clusters(cds_10X)
sce_10X$monocle_clusters<-cds_10X$Cluster
VlnPlot(sce_10X,features=c("nFeature_RNA","nCount_RNA","mitopercent"),group.by="monocle_clusters")#Cells of Cluster8 and Cluster11 show low gene number and UMI number but high mitochondrial gene percentage, indicating that they are dying cells.
cds_10X <- cds_10X[,!(cds_10X$Cluster=="8"| 
                      cds_10X$Cluster=="11")]
sce_10X <- sce_10X[,!(sce_10X$monocle_clusters=="8"|
                      sce_10X$monocle_clusters=="11")]

#3-8 Clustering (Monocle)
cds_10X <- estimateSizeFactors(cds_10X)
cds_10X <- estimateDispersions(cds_10X)
disp_table <- dispersionTable(cds_10X)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_10X <- setOrderingFilter(cds_10X, unsup_clustering_genes$gene_id)
cds_10X <- reduceDimension(cds_10X, max_components = 2, num_dim = 6,
                           reduction_method = 'tSNE', verbose = T)
cds_10X <- clusterCells(cds_10X, num_clusters = 11)
plot_cell_clusters(cds_10X)

#3-9 Trajectory analysis (Monocle)
diff_test_res <- differentialGeneTest(cds_10X,fullModelFormulaStr = "~Cluster",cores = 10)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds_10X <- setOrderingFilter(cds_10X, ordering_genes)
cds_10X <- reduceDimension(cds_10X, max_components = 2, method='DDRTree')
cds_10X <- orderCells(cds_10X) 
cds_10X <- orderCells(cds_10X,root_state = 2)#State1,2,3 represent nEPC state, bGPC state and nNSC-NB state, respectively.
plot_cell_trajectory(cds_10X)

#3-10 Differential expression analysis on branches (Monocle)
Branch_diff <- BEAM(cds_10X)
Branch_diff<- Branch_diff[order(Branch_diff$qval),]

#3-11 Differential expression analysis on clusters (Seurat)
sce_10X <- NormalizeData(sce_10X) 
sce_10X <- FindVariableFeatures(sce_10X, nfeatures = 2000)
sce_10X <- ScaleData(sce_10X,vars.to.regress ="nCount_RNA") 
sce_10X@active.ident<-sce_10X$monocle_clusters
#3-11-1 Marker gene identification for NSC clusters
sce_10X_NSC<-sce_10X[,sce_10X$monocle_clusters ==3|
                      sce_10X$monocle_clusters ==5|
                      sce_10X$monocle_clusters ==6|
                      sce_10X$monocle_clusters ==8]
NSC_cluster_markers<-FindAllMarkers(sce_10X_NSC,only.pos = TRUE)
#3-11-2 Marker gene identification for EPC clusters
sce_10X_EPC<-sce_10X[,sce_10X$monocle_clusters ==7|
                      sce_10X$monocle_clusters ==9]
EPC_cluster_markers<-FindAllMarkers(sce_10X_EPC,only.pos = TRUE)

#3-12 GO enrichment analysis (clusterProfiler)
nNSC_early_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1 
                                               & NSC_cluster_markers$cluster==8]
nNSC_late_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1
                                              & NSC_cluster_markers$cluster==3]
NB1_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>0.7
                                        & NSC_cluster_markers$cluster==6]
NB2_markers <- NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1
                                        & NSC_cluster_markers$cluster==5]
EPC_early_markers <- EPC_cluster_markers$gene[EPC_cluster_markers$avg_log2FC>1
                                              & EPC_cluster_markers$cluster==7]
EPC_late_markers <- EPC_cluster_markers$gene[EPC_cluster_markers$avg_log2FC>1
                                             & EPC_cluster_markers$cluster==9]
NSC_early_markers_BP <- enrichGO(gene=NSC_early_markers,
                                 OrgDb=org.Mm.eg.db,ont='BP',
                                 keyType = 'SYMBOL')
NSC_late_markers_BP <- enrichGO(gene=NSC_late_markers,
                              OrgDb=org.Mm.eg.db,ont='BP',
                              keyType = 'SYMBOL')
NB1_markers_BP <- enrichGO(gene=NB1_markers,
                         OrgDb=org.Mm.eg.db,ont='BP',
                         keyType= 'SYMBOL')
NB2_markers_BP <- enrichGO(gene=NB2_markers,
                         OrgDb=org.Mm.eg.db,ont='BP',
                         keyType= 'SYMBOL')
EPC_early_markers_BP <- enrichGO(gene=EPC_early_markers,
                               OrgDb=org.Mm.eg.db,ont='BP',
                               keyType = 'SYMBOL')
EPC_late_markers_BP <- enrichGO(gene=EPC_late_markers,
                              OrgDb=org.Mm.eg.db,ont='BP',
                              keyType = 'SYMBOL')

#3-13 Trajectory analysis (destiny)
sce_10X$State <- cds_10X$State
ct <- as.SingleCellExperiment(sce_10X)
dm <- DiffusionMap(ct,n_pcs = 50)
