#t-SNE plot(Fig.1D, Fig.S1C-D, Fig.2E)
library(Seurat)
library(monocle)
library(ggplot2)

DimPlot(sce_10X_integrated, reduction = "tsne.harmony",group.by = "harmony_clusters",pt.size = 0.2)+ 
  scale_color_manual(values=palette(c("#98F7E4","#5F72E2","#985AB1",
                                      "#E6A0EA","#FDC471","#A9E0FF",
                                      "#A8E9B5","#DEF3A0","#B83332",
                                      "#5031A7","#967AE1","#B3B2BA",
                                      "#FF947C","#FF5C90")))
DimPlot(sce_10X_integrated, reduction = "tsne.harmony",group.by = "orig.ident",pt.size = 0.2)
plot_cell_clusters(cds_plate,color_by="orig.ident")+
  scale_color_manual(values=palette(c("#00BF7D","#FF62BB")))
plot_cell_clusters(cds_Plate,color_by="Cluster")+ 
  scale_color_manual(values=palette(c("#A9E0FF","#5F72E2","#967AE1",
                                      "#E6A0EA","#B83332","#B3B2BA",
                                      "#FDC471","#DEF3A0","#A8E9B5","#98F7E4")))
FeaturePlot(sce_10X_integrated[,sce_10X_integrated$harmony_clusters=="6"], reduction = "tsne.harmony", features = "Cspg4",pt.size = 0.6)& 
  scale_color_viridis_c()#using Cspg4 gene as an example
plot_cell_clusters(cds_10X,markers="Cspg4",cell_size=0.6)+               
  facet_wrap(~Cluster, nrow = 3)

#Trajectory plot(Fig.1E-G, Fig.S1E, Fig.S1G, Fig.S2B, Fig.3D, Fig.S4E, Fig.S5A, Fig.S5C, Fig.S6A )
library(monocle)
library(ggplot2)
plot_cell_trajectory(cds_10X_integrated, cell_size=0.8,theta=180,show_branch_points = FALSE)
plot_cell_trajectory(cds_Plate)
plot_cell_trajectory(cds_10X_integrated,color_by="harmony_clusters", cell_size=0.8,theta=180,show_branch_points = FALSE)+
  scale_color_manual(values=palette(c("#98F7E4","#5F72E2","#985AB1",
                                      "#E6A0EA","#FDC471","#A9E0FF",
                                      "#A8E9B5","#DEF3A0","#B83332",
                                      "#5031A7","#967AE1","#B3B2BA",
                                      "#FF947C","#FF5C90")))
plot_cell_trajectory(cds_10X_integrated,markers="Gfap",color_by=NULL,cell_size=0.8,use_color_gradient=TRUE,theta=180,show_branch_points = FALSE) #using Gfap gene as an example

#Trajectory plot(Fig.S1F)
library(destiny)
palette(c("#F8766D","#00BA38","#619CFF"))
plot(dm,1:2,col_by = 'State',grid=NULL)

#Pseudotime plot of marker gene expression(Fig.S1H, Fig.2D, Fig.4A, Fig.4D-E, Fig.5E, Fig.5G, Fig.6A)
library(monocle)
plot_genes_branched_pseudotime(cds_10X_integrated[c("Cspg4","Nkx2-2","Olig1","Olig2","Pdgfra","Sox10")],branch_point=1,ncol=2,cell_size = 0.3)#using Fig.2D as an example

#Pseudotime heatmap(Fig.2A, Fig.3B, Fig.3E, Fig.S3B, Fig.S6C)
library(monocle)
NSC_branch <- (!grepl("1",as.character(cds_10X_integrated$State)))
EPC_branch <- (!grepl("3",as.character(cds_10X_integrated$State)))
plot_pseudotime_heatmap(cds_10X_integrated[aNSC_markers,NSC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[qNSC_markers,NSC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[Motile_ciliopathy_associated_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[CP_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[RS_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[NL_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[DA_genes,EPC_branch],num_clusters = 1)
plot_pseudotime_heatmap(cds_10X_integrated[c("Atp2b4","Stk36","Cfap54","Cfap69","Cfap46","Spef2"),NSC_branch],num_clusters = 1,show_rownames = TRUE)
plot_pseudotime_heatmap(cds_10X_integrated[c("Rsph4a","Ropn1l","Nme5","Rsph9"),NSC_branch],num_clusters = 1,show_rownames = TRUE)
plot_pseudotime_heatmap(cds_10X_integrated[c("Ttc28","Ccdc114","Dnah5","Dnah7b","Wdr78","Dnaic1","Dnaic2","Dynlrb2","Dnal1"),NSC_branch],num_clusters = 1,show_rownames = TRUE)
plot_pseudotime_heatmap(cds_10X_integrated[c("Efcab2","Drc1","Iqcg"),NSC_branch],num_clusters = 1,show_rownames = TRUE)
Lysosome_genes_heatmap<-plot_pseudotime_heatmap(cds_10X_integrated[Lysosome_genes_with_CLEAR_motif,EPC_branch],num_clusters = 3,return_heatmap = TRUE)
geneCluster <- cutree(Lysosome_genes_heatmap[["tree_row"]],k=3)
table(geneCluster)

#Branched pseudotime heatmap(Fig.4B-C, Fig.S4A)
library(monocle)
##EPC-fate specific genes
branch_diff_heatmap <- plot_genes_branched_heatmap(cds_10X_integrated[row.names(Branch_diff[1:3000,])],num_clusters = 3 ,return_heatmap = TRUE,cores = 8)
geneCluster2 <- cutree(branch_diff_heatmap[["ph"]][["tree_row"]],k=3)
table(geneCluster2)
EPC_fate_specific_genes <- names(geneCluster2[geneCluster2 ==3])
EPC_fate_specific_genes_non_ciliary <- EPC_fate_specific_genes[grep(FALSE,EPC_fate_specific_genes %in% Cilium_genes)]
##Hydrocephalus genes
hydrocephalus_genes_non_ciliary <- Hydrocephalus_genes[grep(FALSE,Hydrocephalus_genes %in% Cilium_genes)]
hydrocephalus_genes_non_ciliary_heatmap <- plot_genes_branched_heatmap(cds_10X_integrated[hydrocephalus_genes_non_ciliary],return_heatmap = TRUE,cores = 6)
geneCluster3 <- cutree(hydrocephalus_genes_non_ciliary_heatmap[["ph"]][["tree_row"]],k=3)
EPC_hydrocephalus_genes_non_ciliary <- names(geneCluster3[geneCluster3 ==3])
EPC_hydrocephalus_genes_non_ciliary <- Branch_diff[EPC_hydrocephalus_genes_non_ciliary,]
EPC_hydrocephalus_genes_non_ciliary <- rownames(EPC_hydrocephalus_genes_non_ciliary[EPC_hydrocephalus_genes_non_ciliary$qval<0.0001,])
plot_genes_branched_heatmap(cds_10X_integrated[EPC_hydrocephalus_genes_non_ciliary],num_clusters = 1)
##Transcription factors
TF_branch_diff <- Branch_diff[rownames(Branch_diff) %in% TF_genes,]
TF_branch_diff_heatmap<-plot_genes_branched_heatmap(cds_10X_integrated[row.names(TF_branch_diff[TF_branch_diff$qval<0.0001,])],num_clusters = 3,show_rownames = TRUE,return_heatmap = TRUE,cores = 6)
geneCluster4 <- cutree(TF_branch_diff_heatmap[["ph"]][["tree_row"]],k=3)
table(geneCluster4)
EPC_fate_specific_TFs <- names(geneCluster4[geneCluster4 ==2|geneCluster4 ==3])

#Heatmap of marker gene expression(Fig.2C, Fig.3A)
library(Seurat)
library(dplyr)
##NSC clusters
top20_NSC_cluster_markers <- c(head(NSC_cluster_markers[NSC_cluster_markers$cluster==4,],20)$gene,head(NSC_cluster_markers[NSC_cluster_markers$cluster==6,],20)$gene,head(NSC_cluster_markers[NSC_cluster_markers$cluster==12,],20)$gene,head(NSC_cluster_markers[NSC_cluster_markers$cluster==2,],20)$gene)
sce_10X_NSC <- ScaleData(sce_10X_all_NSC,vars.to.regress="nCount_RNA", features=top20_NSC_cluster_markers)
DoHeatmap(sce_10X_NSC, features = top20_NSC_cluster_markers)
##EPC clusters
top50_EPC_cluster_markers <- c(head(EPC_cluster_markers[EPC_cluster_markers$cluster==8,],50)$gene,head(EPC_cluster_markers[EPC_cluster_markers$cluster==9,],50)$gene)
sce_10X_EPC <- ScaleData(sce_10X_all_EPC,vars.to.regress="nCount_RNA", features=top50_EPC_cluster_markers)
DoHeatmap(sce_10X_EPC, features = top50_EPC_cluster_markers)

#Violin plot(Fig.S1B)
library(Seurat)
library(ggplot2)
sce_combined <- merge(sce_10X_all,y=sce_Plate)
sce_combined$lognCount_RNA <- log10(sce_combined$nCount_RNA)
VlnPlot(sce_combined,features="nFeature_RNA",group.by="orig.ident")
VlnPlot(sce_combined,features="lognCount_RNA",group.by="orig.ident")
VlnPlot(sce_combined,features="mito.percent",group.by ="orig.ident")

#Bar plot of GO enrichment(Fig.S2A, Fig.S3A, Fig.S4B, Fig.S4D, Fig.6H)
#NSC_early_markers_BP is used as an example.
library(ggplot2)
library(stringr)
nNSC_early_markers_BP <- nNSC_early_markers_BP[1:12,]
nNSC_early_markers_BP$Description <- str_trunc(nNSC_early_markers_BP$Description,
                                            width=40,side = "right")
nNSC_early_markers_BP$Description <- factor(nNSC_early_markers_BP$Description,
                                         levels=rev(nNSC_early_markers_BP$Description))
ggplot(nNSC_early_markers_BP) +
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

#Venn diagram(Fig.S4C)
library(VennDiagram)
venn_list <- list(Cilia=Cilium_genes, EPC= EPC_fate_specific_genes)
venn.diagram (venn_list, imagetype = "png", fill = c("blue", "green"), 
              alpha = c(0.5, 0.5), filename = "venn.png")

#PCA plot(Fig.S6D)
library(DESeq2)
plotPCA(ntd, intgroup=c("genotypes"))

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

#Bar plots and Pie plots of immunofluorescence results, Bar plots of immunoblot results and qPCR results, together with dot plots of FPKM, were generated using GraphPAD Prism.