# Unraveling Lineage Roadmaps and Fate Determinants to Postnatal Neural Stem Cells and Ependymal Cells in the Developing Ventricular Zone
## Indroduction
The repository contains the code (in `Scripts` folder) for the data processing and analysis to reproduce the results from the paper titled “Unraveling Lineage Roadmaps and Fate Determinants to Postnatal Neural Stem Cells and Ependymal Cells in the Developing Ventricular Zone”(under review). The reference gene sets cited by this paper can be found in `Reference_genesets` folder and the key gene sets identified in this paper are listed in `Output_genesets` folder. 
## Raw sequencing data
The raw data (fastq files) are deposited into ArrayExpress under the accession number E-MTAB-13855 (ChIP-seq), E-MTAB-13856 (bulk RNA-seq) and E-MTAB-13858 (scRNA-seq).
## Softwares/Packages prerequisite
### 1, General programs
R v4.0.3
### 2, R packages
Seurat(v5.2.1)\
DoubletFinder(v2.0.3)\
tidyverse (v1.3.1)\
patchwork (v1.1.1)\
SingleR (v1.4.1)\
monocle (v2.18.0)\
dplyr (v1.1.4)\
plyr (v1.8.6)\
Venndiagram (v1.6.20)\
clusterProfiler (v3.18.1)\
org.Mm.eg.db (v3.12.0)\
destiny (v3.4.0)\
Biobase (v2.50.0)\
ggplot2 (v3.3.5)\
gridExtra (v2.3)\
DESeq2 (v1.30.0)\
stringr (v1.4.0)
### 3, Other bioinformatics utilities
STAR (v2.7.9a)\
fastp (v0.22.0)\
histat2 (v2.1.0)\
picard (v2.23.8)\
samtools (v1.9)\
macs2 (v2.2.7.1)\
homer (v4.11.1)\
bedtools (v2.29.2)\
deeptools (v3.5.0)\
htseq (v0.13.5)\
featureCount (v2.0.6)\
fetchChromSizes, bedGraphToBigWig from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)
## Reproducing the results from the paper
In `Scripts` folder, `scRNA-seq_10X_pipeline`, `scRNA-seq_Plate_pipeline`, `bulk_RNA-seq_pipeline` and `ChIP-seq_pipeline` are provided to perform data processing starting from fastq files. Then the resulting R objects, together with the provided files in `Reference_genesets`, are subject to `Plot_pipeline` to yield the graphs present in the figures and the files stored in `Output_genesets` folder. 
## Contact
zhengjianqun@ojlab.ac.cn
