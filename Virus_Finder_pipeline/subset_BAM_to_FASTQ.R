suppressMessages(library(Rsamtools))
subset_bam<-function(bam_file,read_names_file, outfile){
 bam = scanBam(bam_file)
 reads = scan(read_names_file,what="char",quiet=T)
 bam <- unname(bam)
 what <- c("qname","seq","qual")
 param <- ScanBamParam(what=what)
 elts <- setNames(bamWhat(param), bamWhat(param))
 lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
 bam_df = do.call("DataFrame", lst)
 bam_df$plus=rep("+",nrow(bam_df))
 bam_df$qname1 = paste0("@",bam_df$qname)
 write.table(bam_df[bam_df$qname %in% reads,c(5,2,4,3)], file = outfile, sep = "\n",quote=F, row.names=F, col.names=F)
}

.unlist <- function (x)
{
 x1 <- x[[1L]]
 if (is.factor(x1)) {
 structure(unlist(x), class = "factor", levels = levels(x1))
 } else {
 do.call(c, x)
}
}
