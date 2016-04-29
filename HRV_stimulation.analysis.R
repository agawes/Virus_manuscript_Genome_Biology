############### This is the R code used for analysis of HRV stimulated ALIs ###############
############### used in virus paper submitted to Genome Biology 			###############

### R version 3.1.0

### load libraries
library("DESeq2")
library("edgeR")

### read Ampliseq count table, remove ERCC genes 
data = read.table("HRV_stim.count_table.txt",h=T, row.names=1)

### downsample to the same number of reads per pair for fair comparison #####

downsample_to = colSums(data)
downsample_to[2] = downsample_to[1]
downsample_to[3] = downsample_to[4]
downsample_to[5] = downsample_to[6]

### create a list with vectors representing read to gene mapping
read2gene = list()
for(i in 1:ncol(data)){
	all_reads = unname(unlist(apply(data.frame(Genes=rownames(data),Sample=data[,i] ), 1,function(x){rep(x[1], x[2])})))
	read2gene<-c(read2gene,list(all_reads))
}
names(read2gene) = colnames(data)

#### downsample data to the same number of reads within a pair
data_sub = data.frame(Gene = rownames(data))
for(i in 1:ncol(data)){
	if (colSums(data)[i] != downsample_to[i]){
		subsample = sample(read2gene[[i]], downsample_to[i], replace=F)
		subsample_table = data.frame(table(subsample))
	} else { 
		subsample_table = data.frame(Gene=rownames(data),data=data[,i])
	}
	names(subsample_table) = c("Gene",colnames(data)[i])
	data_sub = merge(data_sub, subsample_table, all.x=T)
} 
## convert NA's to 0's
data_sub[is.na(data_sub)] <- 0
rownames(data_sub) = data_sub$Gene
data_sub = data_sub[,c(2:ncol(data_sub))]

write.table(data_sub, file="HRV_stim.subsampled.count_table.042916.txt",sep="\t",quote=F)

############### run DE between HRV stimulated samples and controls ############### 
############### control for pairing

design = data.frame(row.names = colnames(data_sub), condition = rep(c("HRV","Control"),3), subject = c("ALI_1","ALI_1","ALI_2","ALI_2","ALI_3","ALI_3"))
dds <- DESeqDataSetFromMatrix(
  countData = data_sub,
  colData = design,
  design = ~ subject + condition)
dds$condition <- relevel( dds$condition, "Control" )
dds <- DESeq(dds)
res = results(dds)

### how many DEGs we detect at FDR <5% ?
length(which(res$padj<0.05)) ### 1857

write.table(res, file = "HRV_ampliseq.DEGs.DESeq2.txt", quote = F, sep = "\t", col.names=NA)

############### save top 50 in vitro genes to use in the heatmap in Figure 4  ###############

in_vitro = data.frame(res)
write.table(rownames(head(in_vitro[order(in_vitro$padj),],n=50)), file="Top50_in_vitro_stim_genes.042916.txt",quote=F,row.names=F,col.names=F)


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
col[intersect(which(in_vivo[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1))] <- rgb(255,0,0,95,maxColorValue=255)

### plot comparing fold changes from both experiments
### mark genes up in vivo, but not in vitro in red

pdf("in_vitro_vs_in_vivo.log2FC.042916.pdf")
plot(in_vitro[common_genes,]$log2FoldChange, in_vivo[common_genes,]$log2FoldChange, col = col,pch=16, xlab="in vitro log2FoldChange", ylab = "in vivo log2FoldChange")
abline(v=c(-1,1),lty=2)
abline(h=c(-1,1),lty=2)
dev.off()


length(intersect(which(in_vivo[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1)))  ### 178
in_vivo_genes = rownames(in_vitro)[intersect(which(in_vivo[common_genes,]$log2FoldChange >= 2), which(in_vitro[common_genes,]$log2FoldChange<1))]
write.table(in_vivo_genes, file="in_vivo_genes.042916.txt", quote=F,col.names=F, row.names=F)

save.image("HRV_stim.042916.Rdata")