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

#-----------------------------------------------------------------------------
#Step1: Run Methimpute for cytosine regions
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/globFun.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/MethimputeReg.R"))

#Input files required:
##Folder containing Cytosine regions as Rdata files 
Regionsfolder <- "/Users/rashmi/basedir/DMRcaller/PRODUCED"
##Output directory
myoutput <- "/Users/rashmi/basedir/DMRcaller/test/"
##Folder containing Methimpute outputs 
Methfiles <- "/Users/rashmi/basedir/DMRcaller/methimpute-out"

#Run (U,M,I) 3-state call with include.intermediate=TRUE
#update = one of c("independent","constrained").

runMethimputeRegions(Regionfiles=Regionsfolder,
                     Methfiles=Methfiles,
                     context=c("CG","CHG","CHH"),
                     include.intermediate=FALSE,
                     probability="constrained",
                     out.dir=myoutput)

#-----------------------------------------------------------------------------
#Step2: Extract patterns
#-----------------------------------------------------------------------------

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/methpatterns.R"))

##Output directory
methout <- "/Users/rashmi/basedir/DMRcaller/test"
methpatterns(methout=methout, 
             context=c("CG","CHG","CHH"),
             chr <- c("chr1","chr2","chr3"),
             #chr <- c("chr1","chr2","chr3","chr4","chr5"),
             out.dir=methout, 
             WT="SRR534177")

#-----------------------------------------------------------------------------
# Merging frequency of patterns by contexts into one dataframe (OPTIONAL)

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
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/methpatterns.R"))

fvals <- fread("/Users/rashmi/basedir/DMRcaller/test/chr1_CG_vals.txt", header = TRUE, sep = "\t")
fdensity <- fread("/Users/rashmi/basedir/DMRcaller/test/chr1_CG_meth-patterns-freq.txt",header = TRUE, sep = "\t")
samplesnames <-fread("/Users/rashmi/basedir/DMRcaller/makeRegScripts/DMRs/sample-names.txt",header = FALSE, sep = "\t")

mydf <- fvals[,4:ncol(fvals)]

#mean of all regions that belong to each pattern
mydf.1 <- aggregate(. ~ Pattern.int, data = mydf, FUN = mean)
#replace sample name to a general name
setnames(mydf.1, old = samplesnames$V1, new = samplesnames$V2)

rownames(mydf.1) <- mydf.1$Pattern.int
data_subset <- mydf.1[,2:(ncol(mydf.1))]
data_subset <- t(apply(data_subset, 1, diffr))
data_subset <- data_subset[,-c(NCOL(data_subset))]

my_hclust <- hclust(dist(data_subset), method = "complete")
#plot dendrogram and inspect how many clusters are optimal
#as.dendrogram(my_hclust_pat) %>% plot(horiz = TRUE)
#plot(my_hclust_pat)
#get_clust <- cutree(tree = my_hclust, k = 5)
#get_clust
breaksList <- seq(from=min(data_subset), to=max(data_subset), length.out = 10)
h <- pheatmap(data_subset)
#ordering dataframe of density of patterns to add it as an annotation row to the heatmap 
myrow <- rownames(data_subset[h$tree_row[["order"]],])
fdensity1 <- fdensity[,c("Pattern.int","CG-density")]
fdensity1 <- data.frame(fdensity1)
mynew <- fdensity1[match(myrow, fdensity1$Pattern.int),]
rownames(mynew) <- mynew[,1]
mynew[,1] <- NULL

pheatmap(data_subset, 
         cluster_cols=FALSE,
         #gaps_col = 1,
         breaks=breaksList,
         cluster_rows = my_hclust,
         clustering_distance_rows = "euclidean", 
         #clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         annotation_row = mynew,
         color = colorRampPalette(brewer.pal(n = 3, name = "YlOrRd"))(length(breaksList)))



