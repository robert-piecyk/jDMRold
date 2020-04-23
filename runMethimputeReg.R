#library(devtools)
#devtools::install_github("jlab-code/jDMR")
rm(list=ls())
library(data.table)
library(dplyr)
library(methimpute)
library(stringr)
library(mgcv)
library(ggplot2)
library(plyr)
library(pheatmap)
library(RColorBrewer)

#Run Methimpute for regions
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/globFun.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/MethimputeReg.R"))

Regionsfolder <- "/Users/rashmi/basedir/DMRcaller/PRODUCED"
myoutput <- "/Users/rashmi/basedir/DMRcaller/test/"
Methfiles <- "/Users/rashmi/basedir/DMRcaller/methimpute-out"
runMethimputeRegions(Regionfiles=Regionsfolder,
                     Methfiles=Methfiles,
                     context=c("CG","CHG","CHH"),
                     out.dir=myoutput)

#Extract patterns
#-----------------------------------------------------------------------------

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/methpatterns.R"))
methout <- "/Users/rashmi/basedir/DMRcaller/test"
methpatterns(methout=methout, 
             context=c("CG","CHG","CHH"),
             chr <- c("chr1","chr2","chr3"),
             #chr <- c("chr1","chr2","chr3","chr4"),"chr5",
             out.dir=methout, 
             WT="SRR534177")

#-----------------------------------------------------------------------------
# Merging frequency of patterns outputs by contexts

chr <- c("chr1","chr2","chr3")
myoutput <- "/Users/rashmi/basedir/DMRcaller/test"
for (i in 1:length(chr)){
  pattern=paste0("^",chr[i],".*\\patterns-freq.txt$")
  myfiles <- list.files(myoutput, pattern=pattern, full.names = TRUE)
  filelist <- lapply(myfiles, function(x){
    file <- fread(x)
    return(file)
  })
  df2 <- Reduce(function(x, y) {
    dplyr::full_join(x, y, by=c("Pattern.int","pattern"))
  }, filelist)
  df2 <- df2 %>% dplyr::select("Pattern.int","pattern",everything())
  #fwrite(x=df2, file=paste0(myoutput, "/", chr[i], "_pattern-freq.txt"), quote=FALSE, 
  #       row.names=FALSE, col.names=TRUE, sep="\t")
}
#-----------------------------------------------------------------------------
#make heatmaps to visualize patterns

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/methpatterns.R"))
f1 <- "/Users/rashmi/basedir/DMRcaller/test/chr1_CG_vals.txt"
file <- fread(f1)
mydf <- file[,4:ncol(file)]
#mean of all regions that belong to each pattern
mydf.1 <- aggregate(. ~ Pattern.int, data = mydf, FUN = mean)
rownames(mydf.1) <- mydf.1$Pattern.int
data_subset <- mydf.1[,2:(ncol(mydf.1))]
data_subset <- t(apply(data_subset, 1, diffr))
data_subset <- data_subset[,-c(NCOL(data_subset))]
#first heatmap for the actual values
my_hclust <- hclust(dist(data_subset), method = "complete")
#plot dendrogram and inspect how many clusters are optimal
#as.dendrogram(my_hclust_pat) %>% plot(horiz = TRUE)
#plot(my_hclust_pat)
get_clust <- cutree(tree = my_hclust, k = 5)
get_clust
breaksList <- seq(from=min(data_subset), to=max(data_subset), length.out = 10)
h <- pheatmap(data_subset, 
         cluster_cols=FALSE,
         #gaps_col = 1,
         breaks=breaksList,
         color = colorRampPalette(brewer.pal(n = 3, name = "YlOrRd"))(length(breaksList)),
         cluster_rows = my_hclust,
         clustering_distance_rows = "euclidean", 
         #clustering_distance_cols = "euclidean", 
         clustering_method = "complete")
h
