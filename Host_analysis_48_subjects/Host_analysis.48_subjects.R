############### This is the R code used for analysis of host transcriptome ###############
############### used in dual RNA-seq of virus and host airway manuscript   ###############
############### submitted to Genome Biology 		                       ###############

# R version: R/3.2.0

library(RDAVIDWebService)
library(org.Hs.eg.db)
library("DESeq2")
library("edgeR")
require(gage)
library(scales)
library("WGCNA")
allowWGCNAThreads()
library(RDAVIDWebService)
library(org.Hs.eg.db)

############# read in the data and metadata  #############

data=read.table("AdditionalFile1.txt",h=T,sep="\t",row.names=1)
colnames(data) = gsub("\\.","-",colnames(data))

meta=read.table("metadata.txt",h=T,sep="\t")
rownames(meta) = meta$Subject_ID
meta=meta[colnames(data),]

asthma = factor(gsub("-.+$","",meta$Subject_ID))
virus = meta$pc_viral_reads
virus[virus>0] <- 1
virus=factor(virus)

#############     MDS plot    #############

design = data.frame(row.names = colnames(data), virus=virus, asthma=asthma)
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = design,
  design = ~asthma+virus)

vsd <- varianceStabilizingTransformation(dds)
vstMat = assay(vsd)

pch=as.numeric(design$asthma)
pch[pch==1] <- 16
pch[pch==2] <- 1

# make size of the symbols reflect the viral depth
cex=rep(1,nrow(meta))
meta$Scaled_viral_depth[is.na(meta$Scaled_viral_depth)] <- 0

for (i in 1:nrow(meta)){
 if (meta$Scaled_viral_depth[i] > 1){
 	cex[i] = cex[i] + log10(meta$Scaled_viral_depth[i])
 }
}

### color by DE groups:

col = as.character(meta$Group)
col[col == "VH"] <- "red"
col[col == "VL"] <- "orange"
col[col == "NV"] <- "blue"
col[is.na(col)] <- "grey"

pdf("Figure4.pdf")
plotMDS(vstMat, col = col, pch=pch, gene.selection="common",cex=cex)
legend(0.5,1.6, legend = c("Virus - high","Virus - low","No virus"), col = c("red","orange","blue"), title="Infection", pch=16, bty="n")
legend(2,1.6, legend = c("Asthma","Control"), pch = c(16,1), title="Asthma",bty="n")
dev.off()

######################      differential expression   ####################################
############         Virus high vs No virus by any method          #######################

design$group = meta$Group

select = c(rownames(design[which(design$group == "VH" | design$group == "NV" ),])) 
data_sub = data[,select]
design_sub = design[select,]
design_sub$group <- droplevels(subset(design_sub$group, design_sub$group != "VL"))

dds <- DESeqDataSetFromMatrix(
  countData = data_sub,
  colData = design_sub,
  design = ~asthma+group)

dds <- DESeq(dds)
res_VH = data.frame(results(dds))
length(which(res_VH$padj<0.05))	## 8126

write.table(res_VH, file = "STable2.txt", quote = F, sep = "\t", col.names=NA)

######################      differential expression   ####################################
#############         Virus low vs No virus by any method          #######################

select = c(rownames(design[which(design$group == "VL" | design$group == "NV" ),])) 
data_sub = data[,select]
design_sub = design[select,]
design_sub$group <- droplevels(subset(design_sub$group, design_sub$group != "VH"))

dds <- DESeqDataSetFromMatrix(
  countData = data_sub,
  colData = design_sub,
  design = ~asthma+group)

dds <- DESeq(dds)
res_VL = data.frame(results(dds))
length(which(res_VL$padj<0.05))	## 100

write.table(res_VL, file = "STable5.txt", quote = F, sep = "\t", col.names=NA)

######################    MDS plot - on 100 VIrus-Low DEGs    ######################

select = c(rownames(design[which(design$group == "VL" | design$group == "NV" ),])) 
pick = which(colnames(vstMat) %in% select)
VL =rownames(res_VL[which(res_VL$padj<0.05),])

pdf("Figure5A.pdf")
plotMDS(vstMat[VL,pick], col = col[pick], pch=pch[pick],cex=cex[pick], main="")
legend(0.5,-0.35, legend = c("Virus - low","No virus"), col = c("orange","blue"), title="Infection", pch=16, bty="n")
legend(1,-0.35, legend = c("Asthma","Control"), pch = c(16,1), title="Asthma",bty="n")
dev.off()


#####################  Immune cell enrichment in 100 Virus-Low DEGs ######################   

###### enrichment of different immune cell types signatures                    ###########
###### data generated at: http://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE3982   ###########
###### by comparison of each cell type to all the others                       ###########

immune_dir = "Immune_cells_signatures/"	
files<-list.files(path=immune_dir,recursive=TRUE)
immune_cells = list()
for(f in 1:length(files)){
	name = gsub(".txt","",files[f])
	df = read.table(paste0(immune_dir,files[f]),header=T)
    immune_cells[[name]] = df
}

### get DEGs into GAGE set format:
VL =rownames(res_VL[which(res_VL$padj<0.05),])
VL_up = rownames(res_VL[which(res_VL$padj<0.05 & res_VL$log2FoldChange >0),])
VL_down = rownames(res_VL[which(res_VL$padj<0.05 & res_VL$log2FoldChange <0),])

geneset_modules<-list()
geneset_modules[["LV_DEGs"]]<-VL
geneset_modules[["LV_DEGs_up"]]<-VL_up
geneset_modules[["LV_DEGs_down"]]<-VL_down

### define gage function - returns just the enrichment p-value
gage_tab<-function(deseq2.res, genesets){
  deseq2.fc=deseq2.res$logFC
  names(deseq2.fc)=deseq2.res$Gene.symbol
  input = deseq2.fc
  input<-input[!is.na(names(input))]
  gage_res <- gage(input, gsets = genesets, ref = NULL, samp = NULL, same.dir = T, set.size = c(10, 6000), full.table=T)
  return(gage_res$greater[,4])
}

modules_immune = data.frame(gage_tab(immune_cells[[1]], geneset_modules))
colnames(modules_immune) = names(immune_cells)[1]
for (i in 2:length(immune_cells)){
	res=data.frame(gage_tab(immune_cells[[i]], geneset_modules))
	modules_immune = merge(modules_immune,res,by="row.names")
	rownames(modules_immune) = modules_immune$Row.names
	modules_immune = modules_immune[,-1]
	colnames(modules_immune)[i] = names(immune_cells)[i]
}

write.table(modules_immune, file="Immune_cells_enrichment.100_LV_DEGs.txt",sep="\t",quote=F)

#### find leading edge genes for the neutrophils enrichment through GSEA: ###### 

rnk_list = immune_cells$Nonactivated_neutrophils[order(immune_cells$Nonactivated_neutrophils$logFC, decreasing=T),c(7,6)]
rnk_list$Gene.symbol = gsub("///.+","", rnk_list$Gene.symbol)
rnk_list = rnk_list[-which(rnk_list$Gene.symbol == ""),]
write.table(rnk_list, file="nonact_neutrophils.rnk",sep="\t",quote=F,col.names=F,row.names=F)
write.table(VL_up, file="VL_up.grp",quote=F,col.names=F,row.names=F)

data.rnk.name="nonact_neutrophils.rnk"
gmt.name="VL_up.grp"
results_path = "GSEA" 
rpt_label = "VL_up.neutrophils"

#run GSEA
command <- paste0("java -cp gsea2-2.2.2.jar -Xmx2014m xtools.gsea.GseaPreranked -rnk ",data.rnk.name,
" -gmx ", gmt.name," -collapse false -scoring_scheme weighted -nperm 1000 -set_max 1000 -set_min 10 -zip_report false",
" -plot_top_x 50 -out ", results_path," -gui false -rpt_label ",rpt_label, " > /dev/null")
system(command)

edge=read.table(Sys.glob(paste0(results_path,"/",rpt_label,".GseaPreranked.*/",gmt.name,".xls")),h=T,sep="\t")
edge = edge[edge$CORE.ENRICHMENT == "Yes",]
nt_genes = as.character(edge$PROBE)

#### make a heatmap of the VL DEGs, mark on the side of the rows gene annotation --> neutrophils and EIF2 signaling
scalebluewhitered <- colorRampPalette(c("royalblue", "blanchedalmond", "darkred"), space = "rgb")(100)
col = rep("grey",nrow(design))
col[which(design$group == "NV") ] <- "royalblue"
col[which(design$group == "VL")] <- "darkorange2"
col[which(design$group == "VH")] <- "darkred"

### genes from eIF2 signaling pathway come from IPA analysis
eif2_pathway = c(VL[grep("^RP",VL)][-2],"EIF2AK2","EIF3D")

col_row = rep("white",length(VL))
col_row[which(VL %in% eif2_pathway)] <- "darkred"
col_row[which(VL %in% nt_genes)] <- "royalblue"

pdf("Figure5B.pdf", height=10,width=10)
heatmap(vstMat[VL,], col=scalebluewhitered, ColSideColors=col, RowSideColors=col_row, scale="row", distfun=function (y) dist(y,method ="manhattan"))
dev.off()

### calculate the correlation of collapsed VL genes with viral depth
group=rep(1,length(VL))
cr = collapseRows(vstMat[VL,], group, VL ,method="ME",connectivityBasedCollapsing=T)
VL_collapsed = as.numeric(cr$datETcollapsed)
cor(VL_collapsed, meta$Scaled_viral_depth, method="spearman") 
# [1] 0.7286694


###########################  WGCNA - on Virus high DEGs #################################

setwd("WGCNA/")

VH =rownames(res_VH[which(res_VH$padj<0.05),])
degData = t(vstMat[VH,])
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

######  find WGCNA soft threshold (power)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
pdf(file = "WGCNA_softThreshold.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

######  based on the plot, we set soft-thresholding power to 7

net = blockwiseModules(degData, power = 7, maxBlockSize = 8500,
                     TOMType = "unsigned", minModuleSize = 50,
                     reassignThreshold = 0, mergeCutHeight = 0.05, detectCutHeight = 0.99,
                     numericLabels = TRUE, pamStage = F, pamRespectsDendro = F,
                     saveTOMs = FALSE, verbose = 3)

table(net$colors)

#    0    1    2    3    4    5    6    7    8    9
# 4748 1615  726  203  189  179  164  117  100   85

######  WGCNA coexpression dendrogram ###### 
pdf("SFigure4.pdf", width = 12, height = 9)
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                  "Module colors",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = labels2colors(net$colors)
datME=moduleEigengenes(degData,mergedColors)$eigengenes

design$group = factor(design$group, levels=c("NV","VL","VH"))
order_v_a = order(design$group, design$asthma)
design = design[order_v_a,]

######  inspect module eigengenes ###### 
col = rep("grey",nrow(design))
col[which(design$group == "NV" & design$asthma == "Control")] <- "cadetblue1"
col[which(design$group == "NV" & design$asthma == "Asthma")] <- "blue"
col[which(design$group == "VL" & design$asthma == "Control")] <- "orange"
col[which(design$group == "VL" & design$asthma == "Asthma")] <- "darkorange2"
col[which(design$group == "VH" & design$asthma == "Control")] <- "red"
col[which(design$group == "VH" & design$asthma == "Asthma")] <- "red4"

rownames(datME) = rownames(degData)

pdf("Module_eigengenes.pdf",height=3,width=7)
par(las=3)
for (i in 1:length(unique(moduleColors))){
	which.module = unique(moduleColors)[i]
	ME=datME[order_v_a, paste("ME",which.module, sep="")]
	barplot(ME, col= col,main=which.module, names.arg=rownames(datME)[order_v_a], cex.names=0.5)
}
dev.off()

######  find hub genes for each module, add gene2module assignments ###### 
######  get connectivity values: ###### 
ADJ1=abs(cor(degData,use="p"))^7
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
Alldegrees1$Module = moduleColors
write.table(Alldegrees1, "STable6.txt",sep="\t",quote=F,col.names=NA)

######  find top 10 hub genes in each module: ###### 
gene2module = data.frame(gene=colnames(degData), module=moduleColors)
top10_conn=list()
for(i in 1:length(unique(moduleColors))){
	g = as.character(gene2module$gene[gene2module$module == unique(moduleColors)[i]])
	df = Alldegrees1[g,]
	df = df[order(df$kWithin, decreasing=T),]
	df = df[1:10,]
	top10_conn[[unique(moduleColors)[i]]] = paste(rownames(df),df$kWithin, sep=":")
}

################# run DAVID functional enrichment on module genes ########################

david<-DAVIDWebService$new(email="wesolowskaa@njhealth.org", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

### initialize data frame to store the results
funcChart = data.frame(matrix(ncol = 13, nrow=1))

modules = unique(gene2module$module)
## start from second module (skip grey)
for (j in 2:length(modules)){ 
	### get a vector of gene names in this module
	genes = as.character(gene2module$gene[gene2module$module == modules[j]])

	### get Entrez gene id's
	x <- org.Hs.egSYMBOL2EG
	mapped_genes <- mappedkeys(x)
	xx <- as.list(x[mapped_genes])
	entrez = rep(0,length(genes))
	for (i in 1:length(genes)){
		if (!is.null(xx[genes[i]][[1]])){entrez[i] = xx[genes[i]][[1]]} 
		else {entrez[i] = "NA"}
	}

	### get the gene list into DAVID
	result<-addList(david, entrez,idType = "ENTREZ_GENE_ID",listName=as.character(modules[j]), listType = "Gene")

	### pick the annotation categories
	setAnnotationCategories(david, c("KEGG_PATHWAY","REACTOME_PATHWAY","BIOCARTA","GOTERM_BP_FAT","GOTERM_CC_FAT","GOTERM_MF_FAT","OMIM_DISEASE"))

	### get term cluster report, only report hits with at least 5 genes included
	func = getFunctionalAnnotationChart(david, count=5)
	if (dim(func)[1] >0) {
		func$module = as.character(modules[j])
		func = func[,c(14,1:5,7:13)]
		### concat results
		names(funcChart) = names(func)
		funcChart = rbind(funcChart,func)
	}
}
write.table(funcChart[-1,], file = "STable7.txt", sep = "\t", quote = F, row.names=F)

########## module enrichment for 70 nasal Th2 genes from Poole et al. 2014  #############

th2_genes = scan("../70_Th2_genes_Poole2014.txt",what="character")

enrichment = rep(1, length(unique(moduleColors)))
for(i in 1:length(unique(moduleColors))){
	A=length(intersect(th2_genes, 	gene2module$gene[gene2module$module==unique(moduleColors)[i]])) ## overlap
	B=length(gene2module$gene[gene2module$module==unique(moduleColors)[i]])## module size
	C=length(gene2module$gene)
	D=length(th2_genes)
	enrichment[i]= phyper(A-1,B,C-B,D,lower.tail=F)
}

names(enrichment) = unique(moduleColors)

#         grey        brown    turquoise          red         pink         blue
# 1.0000000000 0.0003303363 0.9999998295 0.7614931294 0.2128793221 0.2859174848
#      magenta        black       yellow        green
# 1.0000000000 1.0000000000 0.8088043376 1.0000000000

############# module enrichment for genes up only in vivo  ######################
### this is the list of genes upregulated in comparison to the in vitro #########
### HRV stimulation experiment #####

immune_genes = scan("../../HRV_stim_in_vitro/in_vivo_only_genes.txt",what="character")
enrichment = rep(1, length(unique(moduleColors)))
for(i in 1:length(unique(moduleColors))){
	A=length(intersect(immune_genes, 	gene2module$gene[gene2module$module==unique(moduleColors)[i]])) ## overlap
	B=length(gene2module$gene[gene2module$module==unique(moduleColors)[i]])## module size
	C=length(gene2module$gene)
	D=length(immune_genes)
	enrichment[i]= phyper(A-1,B,C-B,D,lower.tail=F)
}
names(enrichment) = unique(moduleColors)
enrichment

#         grey        brown    turquoise          red         pink         blue
# 1.000000e+00 1.000000e+00 1.000000e+00 9.835610e-01 3.663210e-20 1.210533e-57
#      magenta        black       yellow        green
# 1.000000e+00 1.000000e+00 9.247766e-02 1.000000e+00

###### enrichment of different immune cell types signatures in our modules ######
###### data generated at: http://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE3982  #########
###### by comparison of each cell type to all the others  ###########

### get gene2module into GAGE set format:
geneset_modules<-list()
for(module in as.character(unique(gene2module$module))){
  geneset_modules[[module]]<-as.character(gene2module$gene[which(gene2module$module == module)])   
}

### now run gage for all immune cells for all modules:
modules_immune = data.frame(gage_tab(immune_cells[[1]], geneset_modules))
colnames(modules_immune) = names(immune_cells)[1]
for (i in 2:length(immune_cells)){
	res=data.frame(gage_tab(immune_cells[[i]], geneset_modules))
	modules_immune = merge(modules_immune,res,by="row.names")
	rownames(modules_immune) = modules_immune$Row.names
	modules_immune = modules_immune[,-1]
	colnames(modules_immune)[i] = names(immune_cells)[i]
}

write.table(modules_immune, file="STable8.txt",sep="\t",quote=F,col.names=NA)

############# correlation of module eigenegenes with viral depth #################
############# in Virus-high subjects only ##############

select =  which( design$group == "VH" )
viral_depth= meta$Scaled_viral_depth

vir_cor_module = rep(0,ncol(datME))
for (i in 1:length(unique(moduleColors))) {
	vir_cor_module[i] = cor(viral_depth[select], datME[select,i],method="spearman")
}
names(vir_cor_module) = colnames(datME)
vir_cor_module

#    MEblack      MEblue     MEbrown     MEgreen      MEgrey   MEmagenta
#  0.2799015   0.7545170  -0.7374795  -0.4551441   0.7691206   0.3237121
#     MEpink       MEred MEturquoise    MEyellow
#  0.8348366   0.5038227  -0.4308049   0.6961028

 save.image("host_analysis.Rdata")
 