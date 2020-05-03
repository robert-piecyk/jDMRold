rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(tidyr)

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/annotateDMRs.R"))

#supply the folder with DMRs (*state-calls-filtered.txt) or at least the co-ordinates as 3 columns ("seqnames", "start", "end")
wd <- "/Users/rashmi/basedir/DMRcaller/test/DMRs/clean-DMRs"

# annotation files
gff.AT <- "/Users/rashmi/basedir/DMRcaller/Annotations/Arabidopsis_thaliana.TAIR10.43.gff3"
gff.TE <- "/Users/rashmi/basedir/DMRcaller/Annotations/TAIR10_TE.gff3"
gff.pr <- "/Users/rashmi/basedir/DMRcaller/Annotations/TAIR10_promoters.gff3"
#gff.intergenic <- "/Users/rashmi/basedir/DMRcaller/Annotations/intergenic.gff3"

#available annotations
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"

annotateDMRs(gff=c(gff.AT, gff.TE, gff.pr),
             annotation=c("gene","promoters","TE"),
             input=wd)

#------------------------------------------------------------------------------------
# Plot the output of DMR-counts.txt
file <- fread("/Users/rashmi/basedir/DMRcaller/test/DMRs/clean-DMRs/DMR-counts.txt", header=TRUE)
file <- separate(data = file, col = sample, into = c("chr", "context"))
mydf <- reshape2::melt(file, measure.vars=c("all","gene", "promoters", "TE", "multiple.annotations"))
colnames(mydf)[3] <- "Annotation"
ggplot(mydf, aes(fill=Annotation, y=value, x=chr)) + geom_bar(position="stack", stat="identity") + 
  labs(title="DMR annotation", x ="", y = "# DMRs overlapping with annotation") + facet_wrap(vars(context))
#fwrite(x=xx, file=paste0(wd, "/annotations/TargetGenesOverlapa.txt"), quote=FALSE, 
#       row.names=FALSE, col.names=TRUE, sep="\t")
