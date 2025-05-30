#There are two replicates.Replicate1 contains five plates (P0_CD133_RGC_Plate1-5) and replicate2 contains three plates (P0_CD133_RGC_Plate6-8).Here P0_CD133_RGC_Plate1 is used as an example.

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
--readFilesIn P0_CD133_RGC_Plate1_Read2 P0_CD133_RGC_Plate1_Read1 \
--soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 10 \
--soloType CB_UMI_Simple \
--soloCBwhitelist whitelist_plate.csv \
--outSAMattributes CB UB \
--outSAMtype BAM SortedByCoordinate \ 
--soloBarcodeReadLength 0
#The whitellist file can be found in the Mapping folder.
#In the output ¡°Solo.out¡± directory, there are ¡°Gene/raw¡± and ¡°Gene/filtered¡± folders, which contain expression values of each gene in each cell. The files in the ¡°Gene/raw¡± directory contain the gene expression of every barcode in the whitelist. The files in the ¡°Gene/filtered¡± contain similar information with barcodes that have too few reads removed. The files in the ¡°Gene/filtered¡± were imported into R for downstream analysis.

#3 Downstream analysis in R
#3-1 Loading the required R packages
library(Seurat)
library(SingleR)
library(monocle)
library(dplyr)
library(destiny)

#3-2 Filtering endothelial cells, microglia and pericytes/fibroblasts (Seurat and Monocle)
#3-2-1 Creating Seurat object
rep1_plate1_counts <- Read10X(data.dir = '/Solo.out/Gene/filtered')
sce_rep1_plate1 <- CreateSeuratObject(counts=rep1_plate1_counts, 
                                      project = "sce_rep1_plate1")
sce_rep1 <- merge(sce_rep1_plate1,
                y=c(sce_rep1_plate2, sce_rep1_plate3, sce_rep1_plate4, sce_rep1_plate5),
                add.cell.ids = c("plate1", "plate2","plate3","plate4","plate5"), 
                project = "rep1")
#3-2-2 Creating Monocle object
counts <- sce_rep1@assays[["RNA"]]@counts
gene_ann <- data.frame(gene_short_name=row.names(counts),row.names= row.names(counts))
fd <- new("AnnotatedDataFrame",data= gene_ann)
pd <- new("AnnotatedDataFrame",data= sce_rep1@meta.data)
cds_rep1 <- newCellDataSet(counts,phenoData=pd,featureData=fd,
                         expressionFamily=negbinomial.size())
#3-2-3 Normalization
cds_rep1 <- estimateSizeFactors(cds_rep1)
cds_rep1 <- estimateDispersions(cds_rep1)
#3-2-4 Clustering
disp_table <- dispersionTable(cds_rep1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_rep1 <- setOrderingFilter(cds_rep1, unsup_clustering_genes$gene_id)
cds_rep1 <- reduceDimension(cds_rep1,max_components= 2, num_dim = 6,
                          reduction_method = 'tSNE', verbose = T)
cds_rep1 <- clusterCells(cds_rep1, num_clusters = 11)
plot_cell_clusters(cds_rep1)
sce_rep1$monocle_clusters<-cds_rep1$Cluster
#3-2-5 Calculation for mitochondrial gene percentage
MITO <- read.csv("Mouse.Mitocarta3.0.csv",header = FALSE)
counts <- as.matrix(sce_rep1@assays[["RNA"]]@counts)
mitopercent <- apply(counts,2,function(x){sum(x>0&rownames(bb)[x]%in%MITO[,1])}) /sce_rep1$nFeature_RNA*100
sce_rep1$mitopercent<- mitopercent
VlnPlot(sce_rep1,features=c("nFeature_RNA","nCount_RNA","mitopercent"),group.by="monocle_clusters") #No cluster shows abnormal gene number, UMI number or mitochondrial gene percentage.
#3-2-6 Celltype annotation(SingleR)
counts <- GetAssayData(sce_rep1)
mg <- celldex::MouseRNAseqData()
common_genes <- intersect(rownames(counts), rownames(mg)) 
counts <- counts[common_genes,]
mg <- mg[common_genes,]
celltype <- SingleR(test = counts, ref = mg, labels = mg$label.fine) 
sce_rep1$celltype <- celltype$labels
#3-2-7 Filtering
table(sce_rep1$monocle_clusters, sce_rep1$celltype)#Cluster5 is constituted of endothelial cells, microglia and pericytes/fibroblasts, which have been annotated by SingleR.
sce_rep1 <- sce_rep1[,!(sce_rep1$Cluster=="5")]#Cells of Cluster5 are excluded from subsequent analysis.

#3-3 Merging replicate1 and replicate2 data (Seurat)
sce_Plate <- merge(sce_rep1, y = sce_rep2, 
                   add.cell.ids = c("rep1", "rep2"), 
                   project = "sce_Plate ")

#3-4 Clustering (Monocle)
counts<-sce_Plate@assays[["RNA"]]@counts
gene_ann<-data.frame(gene_short_name=row.names(counts),row.names = row.names(counts))
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd <- new("AnnotatedDataFrame",data=sce_Plat@meta.data)
cds_Plate<-newCellDataSet(counts,phenoData=pd,featureData=fd,
                          expressionFamily=negbinomial.size())
cds_Plate <- estimateSizeFactors(cds_Plate)
cds_Plate <- estimateDispersions(cds_Plate)
disp_table <- dispersionTable(cds_Plate)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_Plate <- setOrderingFilter(cds_Plate, unsup_clustering_genes$gene_id)
cds_Plate<-reduceDimension(cds_Plate,max_components=2, num_dim = 6,
                           reduction_method = 'tSNE', verbose = T)
cds_Plate <- clusterCells(cds_Plate, num_clusters = 11)
plot_cell_clusters(cds_Plate)

#3-5 Trajectory analysis (Monocle)
diff_test_res <- differentialGeneTest(cds_Plate,fullModelFormulaStr = "~Cluster",cores = 10)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds_Plate <- setOrderingFilter(cds_Plate, ordering_genes)
cds_Plate <- reduceDimension(cds_Plate, max_components = 2, method='DDRTree')
cds_Plate <- orderCells(cds_Plate) 
plot_cell_trajectory(cds_Plate)

#3-6 Trajectory analysis (destiny)
sce_Plate$State <- cds_Plate$State
ct <- as.SingleCellExperiment(sce_Plate)
dm <- DiffusionMap(ct,n_pcs = 50)