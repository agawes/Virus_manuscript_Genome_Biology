############### This is the R code used for analysis of HRV stimulated ALIs ###############
############### used in virus paper submitted to Genome Biology 			###############

### R version 3.1.0

### load libraries
library("DESeq2")

### read in count table 
data = read.table("AdditionalFile2.txt",h=T, row.names=1)

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
length(which(res$padj<0.05)) ### 1857

write.table(res, file = "STable2.HRV_ampliseq.DEGs.DESeq2.txt", quote = F, sep = "\t", col.names=NA)

############# save top 50 in vitro genes to use in the heatmap in Figure 4  ############
illumina = read.table("../Illumina_20_subjects/SFile1.Illumina_20samples.fpkmMat.txt",h=T,row.names=1)
in_vitro = data.frame(res)
merged = merge(in_vitro,illumina,by="row.names")

write.table(head(merged[order(merged$padj),]$Row.names,n=50), file="Top50.in_vitro_stim_genes.txt",quote=F,row.names=F,col.names=F)


############### comparison to in vivo virus signature ####################

in_vitro = data.frame(res)
### read in DE results from the in vivo study
in_vivo = read.table("../Host_analysis_48_subjects/STable4.Virus_high_vs_No_virus.DEGs.DESeq2.txt",h=T)

### genes detected by both Ampliseq and KAPA with iGenomes
common_genes = intersect(rownames(in_vitro),rownames(in_vivo))
length(common_genes) ## 20546

### color the DEGs coming only from Virus vs Control in green
col = rep(rgb(0,0,0,50,maxColorValue=255), length(common_genes))

### color genes upregulated in vivo, but not in vitro --> red:
col[intersect(which(in_vivo[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1))] <- rgb(255,0,0,95,maxColorValue=255)

### plot comparing fold changes from both experiments
### mark genes up in vivo, but not in vitro in red

pdf("SFigure4.in_vitro_vs_in_vivo.log2FC.pdf")
plot(in_vitro[common_genes,]$log2FoldChange, in_vivo[common_genes,]$log2FoldChange, col = col,pch=16, xlab="in vitro log2FoldChange", ylab = "in vivo log2FoldChange")
abline(v=c(-1,1),lty=2)
abline(h=c(-1,1),lty=2)
dev.off()

length(intersect(which(in_vivo[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1)))  ### 199
df = in_vivo[common_genes,]
in_vivo_genes = rownames(df)[intersect(which(in_vivo[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1))]
write.table(in_vivo_genes, file="in_vivo_only_genes.txt", quote=F,col.names=F, row.names=F)

save.image("HRV_stim.Rdata")
