############### This is the R code used for analysis of HRV stimulated ALIs ###############
############### used in virus paper submitted to Genome Biology 			###############

### R version 3.1.0

### load libraries
library("DESeq2")
library("edgeR")

### read Ampliseq count table, remove ERCC genes 
data = read.table("HRV_stim.count_table.txt",h=T, row.names=1)
data = data[-grep("^ERCC-",rownames(data)),]

############### run DE between HRV stimulated samples and controls ############### 
############### control for pairing

design = data.frame(row.names = colnames(data), condition = rep(c("HRV","Control"),3), subject = c("ALI_1","ALI_1","ALI_2","ALI_2","ALI_3","ALI_3"))
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = design,
  design = ~ subject + condition)
dds$condition <- relevel( dds$condition, "Control" )
dds <- DESeq(dds)
res = results(dds)

### how many DEGs we detect at FDR <5% ?
length(which(res$padj<0.05)) ### 2062

write.table(res, file = "HRV_ampliseq.DEGs.DESeq2.txt", quote = F, sep = "\t", col.names=NA)

############### save top 50 in vitro genes to use in the heatmap in Figure 4  ###############

### I actually save top 51 genes, because one of the gene names from Ampliseq was not 
### found in the Illumina nasal brushin RNA-seq data (generated with iGenomes GTF)

write.table(rownames(head(in_vitro[order(in_vitro$padj),],n=51)), file="Top50_in_vitro_stim_genes.042216.txt",quote=F,row.names=F,col.names=F)


############### comparison to in vivo virus signature ####################

in_vitro = data.frame(res)
### read in DE results from the in vivo study
in_vivo = read.table("Virus.asthma_adj.DEGs.DESeq2.txt",h=T)

### genes detected by both Ampliseq and KAPA with iGenomes
common_genes = intersect(rownames(in_vitro),rownames(in_vivo))
length(common_genes) ## 20546

### color the DEGs coming only from Virus vs Control in green
col = rep(rgb(0,0,0,50,maxColorValue=255), length(common_genes))

### color genes upregulated in vivo, but not in vitro --> red:
col[intersect(which(virus[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1))] <- rgb(255,0,0,95,maxColorValue=255)

### plot comparing fold changes from both experiments
### mark genes up in vivo, but not in vitro in red

pdf("in_vitro_vs_in_vivo.log2FC.042216.pdf")
plot(in_vitro[common_genes,]$log2FoldChange, in_vivo[common_genes,]$log2FoldChange, col = col,pch=16, xlab="in vitro log2FoldChange", ylab = "in vivo log2FoldChange")
abline(v=c(-1,1),lty=2)
abline(h=c(-1,1),lty=2)
dev.off()
