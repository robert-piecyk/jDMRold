rm(list=ls())
library(data.table)
library(plyr)
library(dplyr)
library(mgcv)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#-----------------------------------------------------------------------------
# Make heatmaps to visualize patterns
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/methpatterns.R"))

# rc.meth.lvl data
fvals <- fread("/Users/rashmi/basedir/DMRcaller/test/vals_filtered_0.0005-density.txt", 
               header = TRUE, sep = "\t")

# density/frequency output
fdensity <- fread("/Users/rashmi/basedir/DMRcaller/test/StroudMutants/CG_All_methpatterns-freq.txt",
                  header = TRUE, sep = "\t")
#sample names
samplesnames <-fread("/Users/rashmi/basedir/DMRcaller/test/StroudMutants/Stroud-sample-names.txt",
                     header = FALSE, sep = "\t")
#output directory
output.dir <-"/Users/rashmi/basedir/DMRcaller/test"

#-----------------------------------------------------------------------------
mydf <- fvals[,4:ncol(fvals)]

#mean of all regions that belong to each pattern
mydf.1 <- aggregate(. ~ Pattern, data = mydf, FUN = mean)

#replace sample name to a general name
setnames(mydf.1, old = samplesnames$V1, new = samplesnames$V2)

rownames(mydf.1) <- mydf.1$Pattern
data_subset <- mydf.1[,2:(ncol(mydf.1))]
data_subset <- t(apply(data_subset, 1, diffr))
data_subset <- data_subset[,-c(NCOL(data_subset))]

my_hclust <- hclust(dist(data_subset), method = "complete")
#plot dendrogram and inspect how many clusters are optimal
#as.dendrogram(my_hclust_pat) %>% plot(horiz = TRUE)
#plot(my_hclust_pat)
#get_clust <- cutree(tree = my_hclust, k = 5)
#get_clust
#breaksList <- seq(from=min(data_subset), to=max(data_subset), length.out = NROW(data_subset))

quantile_breaks <- function(xs, n ) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
}

breaksList <- quantile_breaks(data_subset, n = 10)

h0 <- pheatmap(data_subset, show_rownames = FALSE)
#ordering dataframe of density of patterns to add it as an annotation row to the heatmap 
myrow <- rownames(data_subset[h0$tree_row[["order"]],])
fdensity1 <- fdensity[,c("Pattern","density")]
fdensity1 <- data.frame(fdensity1)
mynew <- fdensity1[match(myrow, fdensity1$Pattern),]
rownames(mynew) <- mynew[,1]
mynew[,1] <- NULL

pdf(paste(output.dir, "/heatmap_0.0005-density.pdf", sep=""), width = 16, height = 8)
h <- pheatmap(data_subset, 
              cluster_cols=TRUE,
              #gaps_col = 1,
              breaks=breaksList,
              cluster_rows = my_hclust,
              clustering_distance_rows = "euclidean", 
              #clustering_distance_cols = "euclidean", 
              clustering_method = "complete",
              annotation_row = mynew,
              show_rownames = FALSE,
              fontsize_col = 8,
              color = colorRampPalette(brewer.pal(n = 3, name = "YlOrRd"))(length(breaksList)))
print (h)
dev.off()
#--------------------------------------------------------------
#binary matrix
setnames(fmat, old = samplesnames$V1, new = samplesnames$V2)
m <- ifelse(fmat[,5:(ncol(fmat)-1)] == -1, yes = 0, no = 1)
mm <- uniquecombs(m)
mx <- apply(mm[,c(1:ncol(mm))], 1, paste, collapse="")
mxVec <- mx[attr(mm, "index")] 
datas1 <- cbind(mm, mx)
datas1 <- data.frame(datas1)
colnames(datas1)[ncol(datas1)] <- "Pattern"
data_subset1 <- datas1[,-c(NCOL(datas1))]
data_subset1 <- apply(data_subset1, 2, as.numeric)
rownames(data_subset1) <- datas1$Pattern
data_subset1 <- data_subset1[order(match(rownames(data_subset1), myrow)),]
breaksList1 <- seq(from=min(data_subset1), to=max(data_subset1), length.out = 2)
pdf(paste(output.dir, "/chr1_CG_binary.pdf", sep=""), 
    #colormodel = 'cmyk', 
    width = 8, 
    height = 14)
h1 <- pheatmap(data_subset1, 
               cluster_cols=FALSE,
               cluster_rows=FALSE,
               show_rownames=FALSE,
               legend_breaks = c(0,1),
               color = colorRampPalette(brewer.pal(n = 3, name = "YlOrRd"))(length(breaksList1)))

print(h1)
dev.off()

