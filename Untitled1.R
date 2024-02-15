#There are two fastq files:VZ_P0_10X and VZ_P5_10X. Here VZ_P0_10X file was used as an example.

#1,Building the genome index of mm10
STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir STAR_index/mm10/ \
--genomeFastaFiles mm10.fa \
--sjdbGTFfile gencode.vM10.annotation.gtf
#The gtf file can be found in the mapping folder.

#2,Genome mapping
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
