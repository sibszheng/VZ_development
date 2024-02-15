library(monocle)
library(Seurat)
library(SingleR)
library(dplyr)
library(VennDiagram)
library(clusterProfiler)
library(org.Mm.eg.db)
library(destiny)
library(Biobase)
library(ggplot2)
library(gridExtra)

#Monocle
diff_test_res3 <- differentialGeneTest(VZ_monocle3,fullModelFormulaStr = "~Cluster",cores = 10)
ordering_genes3 <- row.names(subset(diff_test_res3, qval < 0.01))
VZ_monocle3 <- setOrderingFilter(VZ_monocle3, ordering_genes3) 
plot_ordering_genes(VZ_monocle3)
VZ_monocle3 <- reduceDimension(VZ_monocle3, max_components = 2,method='DDRTree')
VZ_monocle3 <- orderCells(VZ_monocle3) 
VZ_monocle3 <- orderCells(VZ_monocle3,root_state = 2) 
plot_cell_trajectory(VZ_monocle3) 
plot_cell_trajectory(VZ_monocle3,color_by = "Cluster",cell_size = 1.0)+scale_color_manual(values = palette(c("#A9E0FF","#5F72E2","#967AE1","#E6A0EA","#B83332","#B3B2BA","#FDC471","#DEF3A0","#A8E9B5","#98F7E4")))
plot_cell_trajectory(VZ_monocle3, markers = "gene name",color_by = NULL,use_color_gradient=TRUE)
plot_cell_clusters(VZ_monocle3,color_by = "Cluster",cell_size = 0.6)+scale_color_manual(values = palette(c("#A9E0FF","#5F72E2","#967AE1","#E6A0EA","#B83332","#B3B2BA","#FDC471","#DEF3A0","#A8E9B5","#98F7E4")))
plot_cell_trajectory(VZ_monocle3, color_by = "Pseudotime")
plot_genes_branched_pseudotime(VZ_monocle3[c("gene name")])
NSC<-(!grepl("1",as.character(VZ_monocle3$State)))
EPC<-(!grepl("3",as.character(VZ_monocle3$State)))
plot_pseudotime_heatmap(VZ_monocle3[aNSC_markers,NSC],num_clusters = 1)
plot_pseudotime_heatmap(VZ_monocle3[qNSC_markers,NSC],num_clusters = 1)
plot_pseudotime_heatmap(VZ_monocle3[Motile_ciliopathy_genes,NSC],num_clusters = 1)
plot_pseudotime_heatmap(VZ_monocle3[CP_genes,EPC],num_clusters = 1)
plot_pseudotime_heatmap(VZ_monocle3[RS_genes,EPC],num_clusters = 1)
plot_pseudotime_heatmap(VZ_monocle3[NL_genes,EPC],num_clusters = 1)
plot_pseudotime_heatmap(VZ_monocle3[DA_genes,EPC],num_clusters = 1)
BEAM_res <- BEAM(VZ_monocle3)
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
branch_diff_heatmap<-plot_genes_branched_heatmap(VZ_monocle3[row.names(BEAM_res[1:3000,])],num_clusters = 3,return_heatmap = TRUE)
geneCluster_branch<- cutree(branch_diff_heatmap[["ph"]][["tree_row"]],k=3)
EPC_diff_10X<-names(geneCluster_branch[geneCluster_branch==1])
venn_list <- list(Cilia=Cilium_genes, EPC=EPC_diff_10X)
venn.diagram (venn_list, imagetype = "png", fill = c("blue", "green"), alpha = c(0.5, 0.5), filename = "venn.png")
EPC_diff_10X_nonciliary<-EPC_diff_10X[grep(FALSE,EPC_diff_10X %in% Cilium_genes)]
TF_BEAM_res<- BEAM_res[rownames(BEAM_res) %in% TF_genes,]
branch_TF_diff_heatmap<-plot_genes_branched_heatmap(VZ_monocle3[row.names(TF_BEAM_res[TF_BEAM_res$qval<0.001,])],num_clusters = 3,return_heatmap = TRUE,show_rownames = TRUE)
plot_pseudotime_heatmap(VZ_monocle3[Lysosome_genes_with_CLEAR_motif[!Lysosome_genes_with_CLEAR_motif %in% c("Cpvl","Ctsk","Ctss","Hpse","Mpo","Neu4")],EPC],num_clusters = 3,return_heatmap = TRUE)
hydrocephalus_genes_non_ciliary<-Hydrocephalus_genes[grep(FALSE,hydrocephalus_genes_all %in% Cilium_genes)]
hydrocephalus_genes_non_ciliary_heatmap<-plot_genes_branched_heatmap(VZ_monocle3[hydrocephalus_genes_non_ciliary],num_clusters = 3,return_heatmap = TRUE)
geneCluster_hydrocephalus_genes_non_ciliary<- cutree(hydrocephalus_genes_non_ciliary_heatmap[["ph"]][["tree_row"]],k=3)
EPC_hydrocephalus_genes_non_ciliary<-names(geneCluster_hydrocephalus_genes_non_ciliary[geneCluster_hydrocephalus_genes_non_ciliary==3])
hydrocephalus_genes_screen<-BEAM_res[EPC_hydrocephalus_genes_non_ciliary,]
hydrocephalus_genes_screen<-hydrocephalus_genes_screen[hydrocephalus_genes_screen$qval<0.001,]
hydrocephalus_genes_screen_heatmap<-plot_genes_branched_heatmap(VZ_monocle3[rownames(hydrocephalus_genes_screen)],num_clusters = 1)

#Seurat
sec.combined3 <- NormalizeData(sec.combined3) 
sec.combined3 <- FindVariableFeatures(sec.combined3, nfeatures = 2000)
sec.combined3 <- ScaleData(sec.combined3,vars.to.regress ="nCount_RNA") 
sec.combined3 <- RunPCA(sec.combined3,ndims.print = 1:5, nfeatures.print = 5)
PCAPlot(sec.combined3,group.by="State")
sec.combined3 <- JackStraw(sec.combined3, num.replicate = 100)
sec.combined3 <- ScoreJackStraw(sec.combined3, dims = 1:20)
JackStrawPlot(sec.combined3, dims = 1:20)
ElbowPlot(sec.combined3)
sec.combined3 <- FindNeighbors(sec.combined3, dims = 1:20)
sec.combined3 <- FindClusters(sec.combined3, resolution = 0.4)
sec.combined3 <- RunTSNE(object = sec.combined3, dims.use = 1:20)
sec.combined3 <- RunUMAP(sec.combined3, reduction = "pca", dims = 1:20)
sec.combined3@active.ident<-sec.combined3$monocle_clusters_filter5
sec.combined3_NSC<-sec.combined3[,sec.combined3$monocle_clusters_filter5==3|sec.combined3$monocle_clusters_filter5==5|sec.combined3$monocle_clusters_filter5==6|sec.combined3$monocle_clusters_filter5==8]
NSC_cluster_markers<-FindAllMarkers(sec.combined3_NSC,only.pos = TRUE)
top20_NSC_cluster_markers<-NSC_cluster_markers %>% group_by(cluster) %>% top_n(n = 20, wt =avg_log2FC)
sec.combined3_NSC <- ScaleData(sec.combined3_NSC,vars.to.regress ="nCount_RNA",features = top20_NSC_cluster_markers$gene)
NSC_cluster_names<-c("mid","late2","late1","early")
names(NSC_cluster_names)<-levels(sec.combined3_NSC)
sec.combined3_NSC<-RenameIdents(sec.combined3_NSC,NSC_cluster_names)
top20_NSC_cluster_markers<-rbind(top20_NSC_cluster_markers[top20_NSC_cluster_markers$cluster==8,],top20_NSC_cluster_markers[top20_NSC_cluster_markers$cluster==3,],top20_NSC_cluster_markers[top20_NSC_cluster_markers$cluster==6,],top20_NSC_cluster_markers[top20_NSC_cluster_markers$cluster==5,])
DoHeatmap(sec.combined3_NSC, features = top20_NSC_cluster_markers$gene) 
sec.combined3_EPC<-sec.combined3[,sec.combined3$monocle_clusters_filter5==7|sec.combined3$monocle_clusters_filter5==9]
EPC_cluster_markers<-FindAllMarkers(sec.combined3_EPC,only.pos = TRUE)
top50_EPC_cluster_markers<-EPC_cluster_markers %>% group_by(cluster) %>% top_n(n = 50, wt =avg_log2FC)
sec.combined3_EPC <- ScaleData(sec.combined3_EPC,vars.to.regress ="nCount_RNA",features = top50_EPC_cluster_markers$gene)
DoHeatmap(sec.combined3_EPC, features = top50_EPC_cluster_markers$gene)
DotPlot(VZ_sce_filter5, features =c("Aldh1l1","Aqp4","Fabp7","Gfap","Slc1a3","Vcam1","Ccdc67","Ccno","Dynlrb2","Hydin","Ift20","Rsph1","Ascl1","Egfr","Sox11","Dcx","Dlx1","Tubb3"),group.by ="State") + coord_flip()
#GO analysis
table(NSC_cluster_markers$cluster[NSC_cluster_markers$avg_log2FC>1])
NSC_early_markers<-NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1 & NSC_cluster_markers$cluster==8]
NSC_mid_markers<-NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1 & NSC_cluster_markers$cluster==3]
NSC_late1_markers<-NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>0.7 & NSC_cluster_markers$cluster==6]
NSC_late2_markers<-NSC_cluster_markers$gene[NSC_cluster_markers$avg_log2FC>1 & NSC_cluster_markers$cluster==5]
ego_NSC_early_markers_BP <- enrichGO(gene = NSC_early_markers, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
ego_NSC_mid_markers_BP <- enrichGO(gene = NSC_mid_markers, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
ego_NSC_late1_markers_BP <- enrichGO(gene = NSC_late1_markers, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
ego_NSC_late2_markers_BP <- enrichGO(gene = NSC_late2_markers, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
barplot(ego_NSC_early_markers_BP,showCategory=10,drop=T)
barplot(ego_NSC_mid_markers_BP,showCategory=10,drop=T)
barplot(ego_NSC_late1_markers_BP,showCategory=10,drop=T)
barplot(ego_NSC_late2_markers_BP,showCategory=10,drop=T)
table(EPC_cluster_markers$cluster[EPC_cluster_markers$avg_log2FC>1])
EPC_early_markers<-EPC_cluster_markers$gene[EPC_cluster_markers$avg_log2FC>1 & EPC_cluster_markers$cluster==7]
EPC_late_markers<-EPC_cluster_markers$gene[EPC_cluster_markers$avg_log2FC>1 & EPC_cluster_markers$cluster==9]
ego_EPC_early_markers_BP <- enrichGO(gene = EPC_early_markers, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
ego_EPC_late_markers_BP <- enrichGO(gene = EPC_late_markers, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
barplot(ego_EPC_early_markers_BP,showCategory=10,drop=T)
barplot(ego_EPC_late_markers_BP,showCategory=10,drop=T)
ego_EPC_diff_10X_BP<-enrichGO(gene = EPC_diff_10X, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
ego_EPC_diff_10X_nonciliary_BP<-enrichGO(gene = EPC_diff_10X_nonciliary, OrgDb = org.Mm.eg.db, ont='BP',keyType = 'SYMBOL')
barplot(ego_EPC_diff_10X_BP,showCategory=10,drop=T)
barplot(ego_EPC_diff_10X_nonciliary_BP,showCategory=10,drop=T)
#Diffusion map
ct <- as.SingleCellExperiment(sce.combined3)
dm <- DiffusionMap(ct,n_pcs = 50)
qplot(y = eigenvalues(dm)) + theme_minimal() +labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')
palette(c("#F8766D","#00BA38","#619CFF"))
plot(dm,1:2,col_by = 'State')
palette(c("#F8766D","#DA950D","#A3A500","#39B600","#02BF7E","#00BFC4","#0AB3F6","#9590FF","#E76BF3","#FF62BC"))
plot(dm,1:2,col_by = 'monocle_clusters')
dpt <- DPT(dm)
grid.arrange(plot(dpt))
dpt <- DPT(dm,tips = c(6151,9782,10661))
grid.arrange(plot(dpt,col_by = 'gene name'))

#EPC stateËæ»úÈ¡Êý
sce.combined3_EPC_state<-sce.combined3[,sce.combined3$State==1]
bb<-as.matrix(sce.combined3_EPC_state@assays[["RNA"]]@counts)
EPC_state_genes<-apply(bb,1,sum)
EPC_state_genes<-EPC_state_genes[EPC_state_genes>1]                   
EPC_state_gene_names<-names(EPC_state_genes)
sampled<-sample(1:22088,1320)
EPC_state_gene_names_sampled<-EPC_state_gene_names[sampled]
table(EPC_state_gene_names_sampled %in% Cilium_genes)
