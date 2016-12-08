suppressMessages(library(ggplot2))
cov_plot<-function(cov_file,blast_contigs,gene_structure,subject,virus){
 cov = read.table(cov_file,h=F,sep="\t")
 names(cov) = c("name","st","e","pos","dp")
 blast = read.table(blast_contigs,h=F)
 bl_vir = blast[grep(virus,blast$V2),]
 bl_vir$V13 = as.numeric(gsub(".+cov_","",bl_vir$V1))
 names(bl_vir) = c("contig","vir","pc","l","r","m","x1","l1","start","end","p","mm","cov")
 genes = read.table(gene_structure,sep="\t")
 names(genes) = c("Start","End","Gene")
 ymin=-0.25*(max(cov$dp)+1)
 ymax=max(cov$dp)+1

 p1 <- ggplot(environment = environment()) + geom_area(data=cov, aes(x=pos, y=dp), fill = "red", alpha =0.5) + ylim(-0.25*(max(cov$dp)+1),max(cov$dp)+1) + theme_bw() + xlab("Position along the genome") + ylab("Depth") + ggtitle(virus)
 p2 <- p1 + geom_segment(data=bl_vir, mapping=aes(x=start, y = cov, xend = end, yend = cov))

 p3 <- p2 + geom_rect(data=genes, mapping=aes(xmin=Start, xmax=End,ymin=-0.25*(max(cov$dp)+1),ymax=-0.025*(max(cov$dp)+1),fill=Gene), colour="black", alpha=0.5) + geom_text(data=genes, aes(x=(Start+End)/2,y=0.6*(-0.25*(max(cov$dp)+1)), label=Gene), size=4) + theme(legend.position="none")
 
 pdf(paste0(subject,"/",subject,".",virus,".pdf"), width=12,height=4)
 plot(p3)
 dev.off()
}
