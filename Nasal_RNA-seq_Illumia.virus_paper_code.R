######## This is the R code used for host analysis of Nasal RNA-seq Illumina data ########
######## used in virus paper submitted to Genome Biology 			###############

### R version 3.0.1

### load libraries

library(cummeRbund)
library(gplots)
library(RColorBrewer)
library(DESeq2)

############### Figure 4 - heatmap of top 50 in vitro HRV stimulated genes ###############

### load cummeRbund object with Cufflinks DE results

cuff_data <- readCufflinks()
fpkmMat <- repFpkmMatrix(genes(cuff_data));

rhino_genes = scan("Top50_in_vitro_stim_genes.042216.txt",what="character")
rhino_fpkms = fpkmMat[as.character(rhino_genes),]
rhino_fpkms = rhino_fpkms[-36,] ### this gene is not present or has a different name in iGenomes GTF
colnames(rhino_fpkms) = c("asthma_0","asthma_1","asthma_2","asthma_3","asthma_4","asthma_5","asthma_6","asthma_7","asthma_8","asthma_9","control_0","control_1","control_2","control_3","control_4","control_5","control_6","control_7","control_8","control_9")

scalegreenblackred <- colorRampPalette(c("darkgreen", "blanchedalmond", "darkred"), space = "rgb")(100)
pdf("HRV_stim_top_genes.heatmap.042216.pdf")
heatmap(data.matrix(log10(rhino_fpkms+1)), col=scalegreenblackred, scale="row")
dev.off()

############### Differential expression - virus vs no-virus #####################
### this was run on count tables generated with HtSeq on the Tophat mapped BAM files
### R version 3.1.0

data = read.table("Nasal_SF.merged_counts.txt",h=T, row.names=1)

virus=rep("no_virus",20)
virus[17] = "virus"
virus[20] = "virus"

design = data.frame(row.names = colnames(data), condition = virus)

dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = design,
  design = ~ condition)

dds$condition <- relevel( dds$condition, "no_virus" )
dds <- DESeq(dds)
res = results(dds)

### find biomarkers of viral infection - genes with the highest log2FoldChanges
### in virus samples compared to non-viral samples
head(res[order(res$log2FoldChange, decreasing=T),])
