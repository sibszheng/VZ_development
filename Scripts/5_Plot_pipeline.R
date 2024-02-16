#t-SNE plot(Fig.1B, Fig.S1C, Fig.2E)
library(monocle)
library(ggplot2)
plot_cell_clusters(cds_10X,color_by="Cluster",cell_size=0.6)+ 
  scale_color_manual(values=palette(c("#A9E0FF","#5F72E2","#967AE1","#E6A0EA","#B83332","#B3B2BA","#FDC471","#DEF3A0","#A8E9B5","#98F7E4")))
plot_cell_clusters(cds_Plate,color_by="Cluster")+ 
  scale_color_manual(values=palette(c("#A9E0FF","#5F72E2","#967AE1","#E6A0EA","#B83332","#B3B2BA","#FDC471","#DEF3A0","#A8E9B5","#98F7E4")))
plot_cell_clusters(cds_10X,markers="Cspg4",cell_size=0.6)+               
  facet_wrap(~Cluster, nrow = 3)#using Cspg4 gene as an example

#Trajectory plot(Fig.1C-E, Fig.S1C, Fig.2D, Fig3D, Fig.5C, Fig.5F, Fig.6A, Fig.S1E, Fig.S4E)
library(monocle)
library(ggplot2)
plot_cell_trajectory(cds_10X, cell_size=1.0)
plot_cell_trajectory(cds_Plate)
plot_cell_trajectory(cds_10X,color_by="Cluster",cell_size=1.0)+
  scale_color_manual(values=palette(c("#A9E0FF","#5F72E2","#967AE1","#E6A0EA","#B83332","#B3B2BA","#FDC471","#DEF3A0","#A8E9B5","#98F7E4")))
plot_cell_trajectory(cds_10X,markers="Gfap",color_by=NULL,cell_size=1.0,use_color_gradient=TRUE) #using Gfap gene as an example

#Trajectory plot(Fig.S1D)
library(destiny)
palette(c("#F8766D","#00BA38","#619CFF"))
plot(dm,1:2,col_by = 'State')

#Pseudotime plot of marker gene expression(Fig.S1F, Fig.4A, Fig.4D-E, Fig.S5A, Fig.S5C, Fig.S6A)
library(monocle)
plot_genes_branched_pseudotime(cds_10X[c("Ccno","Cfap43")],branch_point=1,ncol=2,cell_size = 0.3)#using Ccno gene and Cfap43 gene as examples

#Pseudotime plot of averaged expression of marker genes(Fig.1F)
library(ggplot2)
library(monocle)
library(plyr)
library(dplyr)
marker_genes = list(
  GPCs = c("Gfap","Slc1a3","Vcam1","Aldh1l1","Aqp4","Fabp7"), 
  EPCs = c("Ccdc67","Hydin", "Rsph1", "Dynlrb2", "Ift20", "Tekt4"), 
  NSCs = c("Ascl1", "Egfr", "Sox11"), 
  NBs = c("Dcx","Dlx1", "Tubb3"))
#1 Writing functions
psdt_branch <- function(cds){
  cds_subset <- buildBranchCellDataSet(cds = cds, progenitor_method = 'duplicate')
  CM <- exprs(cds_subset)
  CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(CM)))
  min_expr <- cds_subset@lowerDetectionLimit
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs <- merge(cds_exprs, pData(cds_subset), by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression <- round(cds_exprs$expression)
  new_data<-data.frame(Pseudotime=pData(cds_subset)$Pseudotime,Branch=pData(cds_subset)$Branch)
  full_model_expectation<-genSmoothCurves(cds_subset, cores=1, trend_formula="~ sm.ns(Pseudotime, df=3) * Branch", relative_expr = T, new_data = new_data)
  colnames(full_model_expectation) <- colnames(cds_subset)
  cds_exprs$full_model_expectation<-apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$full_model_expectation[cds_exprs$full_model_expectation< min_expr] <- min_expr
  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  cds_exprs = subset(cds_exprs, select = c("Cell", "expression", "f_id", "Pseudotime", "Branch", "State", "full_model_expectation"))
  return(cds_exprs)}
monocle_theme_opts <- function(){
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())}
psdt_mean_branch <- function(psdt_expr_test){
  psdt_mean_test <- psdt_expr_test %>%
    group_by(Cell) %>%
    reframe(expression = mean(expression),
            Pseudotime = Pseudotime,
            Branch = Branch,
            State = State,
            full_model_expectation = mean(full_model_expectation)) %>%distinct()
  return(psdt_mean_test)}
#2 Geting the expression matrix(GPC is used as an example)
cds_10X_GPC = cds_10X[marker_genes$GPCs,]
cds_10X_GPC_expr = psdt_branch(cds_10X_GPC)
#3 Geting the averaged expression values of the marker genes
cds_10X_GPC_expr_mean = psdt_mean_branch(cds_10X_GPC_expr)
#4 Adding a celltype column
GPCs = data.frame(celltype = rep("GPCs", nrow(cds_10X_GPC_expr_mean)))
cds_10X_GPC_expr_mean = cbind(GPCs, cds_10X_GPC_expr_mean)
#5 combining all data
cds_10X_all_expr_mean=rbind(cds_10X_GPC_expr_mean, cds_10X_EPC_expr_mean,
                            cds_10X_NSC_expr_mean, cds_10X_NB_expr_mean)
ggplot(aes(Pseudotime, expression, color_by = celltype, linetype = Branch),
       data = psdt_mean_sample_B) +
  scale_y_log10() + 
  geom_line(aes(x = Pseudotime, y = full_model_expectation, colour = celltype), 
            data = psdt_mean_sample_B,lwd=1.5) +
  monocle_theme_opts()+ 
  ylab("Expression") + 
  xlab("Pseudotime (stretched)")

#Pseudotime heatmap(Fig.2A, Fig.3B, Fig.3E, Fig.S3B, Fig.S6C)
library(monocle)
NSC_branch <- (!grepl("1",as.character(cds_10X$State)))
EPC_branch <- (!grepl("3",as.character(cds_10X$State)))
plot_pseudotime_heatmap(cds_10X[aNSC_markers,NSC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[qNSC_markers,NSC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[Motile_ciliopathy_associated_genes,NSC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[CP_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[RS_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[NL_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[DA_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X[Lysosome_genes_with_CLEAR_motif,EPC_branch],num_clusters = 3)

#Branched pseudotime heatmap(Fig.4B-C, Fig.S4A)
library(monocle)
##EPC-fate specific genes
branch_diff_heatmap <- plot_genes_branched_heatmap(cds_10X[row.names(Branch_diff[1:3000,])],num_clusters = 3 ,return_heatmap = TRUE)
geneCluster <- cutree(branch_diff_heatmap[["ph"]][["tree_row"]],k=3)
table(geneCluster)
EPC_fate_specific_genes <- names(geneCluster [geneCluster ==1])
EPC_fate_specific_genes_non_ciliary <- EPC_fate_specific_genes[grep(FALSE,EPC_fate_specific_genes %in% Cilium_genes)]
##Hydrocephalus genes
hydrocephalus_genes_non_ciliary <- Hydrocephalus_genes[grep(FALSE,Hydrocephalus_genes %in% Cilium_genes)]
hydrocephalus_genes_non_ciliary_heatmap <- plot_genes_branched_heatmap(cds_10X[hydrocephalus_genes_non_ciliary],num_clusters = 3,return_heatmap = TRUE)
geneCluster <- cutree(hydrocephalus_genes_non_ciliary_heatmap[["ph"]][["tree_row"]],k=3)
EPC_hydrocephalus_genes_non_ciliary <- names(geneCluster[geneCluster ==3])
EPC_hydrocephalus_genes_non_ciliary <- Branch_diff[EPC_hydrocephalus_genes_non_ciliary,]
EPC_hydrocephalus_genes_non_ciliary <- rownames(EPC_hydrocephalus_genes_non_ciliary[EPC_hydrocephalus_genes_non_ciliary$qval<0.001,])
plot_genes_branched_heatmap(cds_10X[EPC_hydrocephalus_genes_non_ciliary],num_clusters = 1)
##Transcription factors
TF_branch_diff <- Branch_diff[rownames(Branch_diff) %in% TF_genes,]
TF_branch_diff_heatmap<-plot_genes_branched_heatmap(cds_10X[row.names(TF_branch_diff[TF_branch_diff$qval<0.001,])],num_clusters = 3,show_rownames = TRUE)
geneCluster <- cutree(TF_branch_diff_heatmap[["ph"]][["tree_row"]],k=3)
table(geneCluster)
EPC_fate_specific_TFs <- names(geneCluster[geneCluster ==2|geneCluster ==3])

#Heatmap of marker gene expression(Fig.2C, Fig.3A)
library(Seurat)
library(dplyr)
##NSC clusters
top20_NSC_cluster_markers <- NSC_cluster_markers %>% group_by(cluster) %>% top_n(n = 20, wt =avg_log2FC)
sce_10X_NSC <- ScaleData(sce_10X_NSC,vars.to.regress="nCount_RNA", features=top20_NSC_cluster_markers$gene)
DoHeatmap(sce_10X_NSC, features = top20_NSC_cluster_markers$gene)
##EPC clusters
top50_EPC_cluster_markers <- EPC_cluster_markers %>% group_by(cluster) %>% top_n(n = 50, wt =avg_log2FC)
sce_10X_EPC <- ScaleData(sce_10X_EPC,vars.to.regress="nCount_RNA", features=top50_EPC_cluster_markers$gene)
DoHeatmap(sce_10X_EPC, features = top50_EPC_cluster_markers$gene)

#Venn diagram(Fig.S4C)
library(VennDiagram)
venn_list <- list(Cilia=Cilium_genes, EPC= EPC_fate_specific_genes)
venn.diagram (venn_list, imagetype = "png", fill = c("blue", "green"), 
              alpha = c(0.5, 0.5), filename = "venn.png")

#Violin plot(Fig.S1B)
library(Seurat)
library(ggplot2)
sce_10X$orig.ident <- NULL
sce_Plate$orig.ident <- NULL
sce_combined <- merge(sce_10X,y=sce_SS3,add.cell.ids=c("10X","Plate"))
sce_combined$lognCount_RNA <- log10(sce_combined$nCount_RNA)
VlnPlot(sce_combined,features="nFeature_RNA",
        group.by="orig.ident",cols=palette(c("#FB8173","#8ED4C8")))+
  scale_y_continuous(limits=c(0,9000),breaks = c(0,2000,4000,6000,8000))+
  theme(legend.position = 'none')
VlnPlot(sce_combined,features="lognCount_RNA",
        group.by="orig.ident",cols = palette(c("#FB8173","#8ED4C8")))+
  scale_y_continuous(limits=c(0,5.4),breaks = c(0,1,2,3,4,5))+
  theme(legend.position = 'none')
VlnPlot(sce_combined,features="mito.percent",
        group.by ="orig.ident",cols = palette(c("#FB8173","#8ED4C8")))+
  scale_y_continuous(limits=c(0,1.4),breaks = c(0,0.5,1))+
  theme(legend.position = 'none')

#Bar plot of GO enrichment(Fig.S2A, Fig.S3A, Fig.6H, Fig.S4B, Fig.S4D)
#NSC_early_markers_BP is used as an example.
library(ggplot2)
library(stringr)
NSC_early_markers_BP <- NSC_early_markers_BP[1:10,]
NSC_early_markers_BP$Description <- str_trunc(NSC_early_markers_BP$Description,
                                            width=40,side = "right")
NSC_early_markers_BP$Description <- factor(NSC_early_markers_BP$Description,
                                         levels=rev(NSC_early_markers_BP$Description))
ggplot(NSC_early_markers_BP) +
  aes(x = Description, y = Count, fill = -log10(p.adjust)) +
  geom_bar(stat = "identity") + 
  scale_fill_continuous(low="#8DC1DE", high="#024395") +
  coord_flip() + 
  theme(axis.text = element_text(size=12,face="plain",color="black"), 
        axis.text.y = element_text(hjust=1,vjust=0.6), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.direction = "vertical", 
        legend.background = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = "black"), 
        plot.background = element_blank())

#Volcano plot(Fig.6F)
library(ggplot2)
df <- as.data.frame(res_TFEB_KD)[,c(2,6)]
df$group <- as.factor(ifelse(df$padj<0.001&abs(df$log2FoldChange)>=0.5, ifelse(df$log2FoldChange>= 0.5 ,'up','down'),'ns'))  
df <- na.omit(df)
ggplot(df, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(color = group),alpha=0.8, size=3)+
  scale_color_manual(values = palette(c("#18499D","#BEBEBE","#DE1F1C")))+
  theme_bw()+
  theme(legend.title = element_blank(),panel.grid = element_blank())

#PCA plot(Fig.S6D)
library(DESeq2)
plotPCA(ntd, intgroup=c("genotypes"))

#Bar plots of immunofluorescence results and qPCR results, together with dot plots of FPKM, were generated using GraphPAD Prism.