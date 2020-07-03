#--------------------------------------------------------
# extract hot- and coldspots in Stroud mutants
#--------------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)
library(GenomicRanges)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

DMean <- fread("/Users/rashmi/basedir/DMRcaller/test/StroudMutants/AB-out/genome-bins-resolutions-3July/MA1_1_divergence-out-10Kb.txt", 
               sep="\t", header=TRUE)
chrArms <- fread("/Users/rashmi/basedir/DMRcaller/test/StroudMutants/myfiles/chr_arms.txt", header=TRUE)
samplenames <- fread("/Users/rashmi/basedir/DMRcaller/test/StroudMutants/myfiles/Stroud-sample-names.txt", 
                     header=FALSE)
CG_status_calls <- fread("/Users/rashmi/basedir/DMRcaller/test/StroudMutants/myfiles/CG_statuscalls.txt")

#replace NA with 0
DMean[is.na(DMean)] <- 0

#--------------------------------------------------------
#extract high and low Dmean values respectively
#--------------------------------------------------------
D.high <- DMean[DMean$meanD > quantile(DMean$meanD, 0.90),]
D.low <- DMean[DMean$meanD < quantile(DMean$meanD, 0.10),]

#create Granges object
D.high.gr <- GRanges(seqnames=D.high$seqnames, 
                     ranges=IRanges(start=D.high$start, 
                                    end=D.high$end),meanD=D.high$meanD)
D.low.gr <- GRanges(seqnames=D.low$seqnames, 
                    ranges=IRanges(start=D.low$start, 
                                   end=D.low$end),meanD=D.low$meanD)
chrArms.gr <- GRanges(seqnames=chrArms$seqnames, 
                      ranges=IRanges(start=chrArms$start, 
                                     end=chrArms$end))
mygr <- GRanges(seqnames=CG_status_calls$seqnames, 
                ranges=IRanges(start=CG_status_calls$start, 
                               end=CG_status_calls$end))

#--------------------------------------------------------
# Overlaps with chromosome arms
#--------------------------------------------------------
overlaps_D.high <- findOverlaps(D.high.gr, chrArms.gr)
bins_D.high <- D.high.gr[queryHits(overlaps_D.high)]

overlaps_D.low <- findOverlaps(D.low.gr, chrArms.gr)
bins_D.low <- D.low.gr[queryHits(overlaps_D.low)]

#--------------------------------------------------------
# Remove overlapping bins between bins_D.high and bins_D.low
#--------------------------------------------------------
unique.bins_D.high <- bins_D.high[!bins_D.high %over% bins_D.low,]
unique.bins_D.low <- bins_D.low[!bins_D.low %over% bins_D.high,]

#-----------------------------------------------------------
# collapse the 400bp bins into hot and cold bins and add ids
#-----------------------------------------------------------
hotspots <- reduce(unique.bins_D.high)
coldspots <- reduce(unique.bins_D.low)
hotspots.id <- data.frame(hotspots) %>% mutate(id = paste0("HS",row_number(),"_chr",seqnames))
coldspots.id <- data.frame(coldspots) %>% mutate(id = paste0("CS",row_number(),"_chr",seqnames))
mcols(hotspots) <- cbind.data.frame(mcols(hotspots), id=hotspots.id$id)
mcols(coldspots) <- cbind.data.frame(mcols(coldspots), id=coldspots.id$id)
hotspots
coldspots

#-------------------------------------------------------------------
# PLOT1: Plot the length distributions of hot- and coldspots
#-------------------------------------------------------------------
df1 <- data.frame(hotspots)
df2 <- data.frame(coldspots)
df1$gp <- "hotspot"
df2$gp <- "coldspot"
fdf <- rbind(df1,df2)

gg <- ggplot(fdf, aes(x=width, fill=gp)) + 
  #geom_histogram(aes(y = ..density..), binwidth=bw, colour="black", 
  #               alpha=0.6, position="identity", size=0.08) +
  geom_density(aes(y = ..density.., colour=gp), 
               alpha = 0.4, position = "identity", size=0.7) +
  scale_x_log10(breaks = scales::log_breaks(n = 12), expand = c(0, 0)) +
  labs(x="log10 (length)")
p <- gg + scale_fill_manual(values=c("red", "blue")) + 
  scale_color_manual(values=c("red","blue")) + 
  theme(strip.background = element_rect(colour = "black", fill="grey"),
        panel.border = element_rect(fill = NA),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) 

print(p)

#--------------------------------------------------------
# Find hot- or cold- bins in Stroud mutants 
#--------------------------------------------------------
overlaps_hotspots_Stroud <- findOverlaps(mygr, hotspots)
hotspots_Stroud <- mygr[queryHits(overlaps_hotspots_Stroud)]
mcols(hotspots_Stroud)$id <- CharacterList(split(hotspots$id[subjectHits(overlaps_hotspots_Stroud)], 
                                                 queryHits(overlaps_hotspots_Stroud)))

overlaps_coldspots_Stroud <- findOverlaps(mygr, coldspots)
coldspots_Stroud <- mygr[queryHits(overlaps_coldspots_Stroud)]
mcols(coldspots_Stroud)$id <- CharacterList(split(coldspots$id[subjectHits(overlaps_coldspots_Stroud)], 
                                                 queryHits(overlaps_coldspots_Stroud)))

#---------------------------------------------------------------------------------
# Extract methylation status calls of hot- or cold- bins in Stroud mutants 
#---------------------------------------------------------------------------------
hotspots_Stroud <- data.frame(hotspots_Stroud)
hot.df <- merge(hotspots_Stroud,CG_status_calls, by=c("seqnames","start","end"))
hot.df <- as.data.frame(lapply(hot.df, unlist))

coldspots_Stroud <- data.frame(coldspots_Stroud)
cold.df <- merge(coldspots_Stroud,CG_status_calls, by=c("seqnames","start","end"))
cold.df <- as.data.frame(lapply(cold.df, unlist))

setnames(hot.df, old=colnames(hot.df[7:NCOL(hot.df)]), new=samplenames$V2)
setnames(cold.df, old=colnames(cold.df[7:NCOL(cold.df)]), new=samplenames$V2)

#------------------------------------------------------------------------------------------------
# Count %methylation status : number of bins with M status over total bins in a hotspot/coldspot
#------------------------------------------------------------------------------------------------
a <- hot.df %>% 
  group_by(id) %>%
  count() %>%
  as.data.frame()

b <- hot.df %>% group_by(id) %>%
  summarise_at(c(6:92), function(x) sum(x == 1)) %>% as.data.frame()

meth.count.hotspots <- merge(a, b, by="id")
perc.meth.count.hotspots <- meth.count.hotspots[,3:NCOL(meth.count.hotspots)] / meth.count.hotspots[,2]

c <- cold.df %>%
  group_by(id) %>%
  count() %>%
  as.data.frame()

d <- cold.df %>% group_by(id) %>%
  summarise_at(c(6:92), function(x) sum(x == 1)) %>% as.data.frame()

meth.count.coldspots <- merge(c,d,by="id")
perc.meth.count.coldspots <- meth.count.coldspots[,3:NCOL(meth.count.coldspots)] / meth.count.coldspots[,2]


#--------------------------------------------------------
# PLOT2: Plot %methylation status of hot- and cold- spots in WT
#--------------------------------------------------------
perc.meth.count.hotspots$gp <- "hotspot"
perc.meth.count.coldspots$gp <- "coldspot"

plot.meth.status <- rbind(perc.meth.count.hotspots[,c("WT","gp")],
              perc.meth.count.coldspots[,c("WT","gp")])

gg <- ggplot(plot.meth.status, aes(x=WT, fill=gp)) + 
  #geom_histogram(aes(y = ..density..), binwidth=bw, colour="black", 
  #               alpha=0.6, position="identity", size=0.08) +
  geom_density(aes(y = ..density.., colour=gp), 
               alpha = 0.4, position = "identity", size=0.7)
#scale_x_log10(breaks = scales::log_breaks(n = 6), expand = c(0, 0)) +
#labs(x="log10 (length)")
p <- gg + scale_fill_manual(values=c("red", "blue")) + 
  scale_color_manual(values=c("red","blue")) + 
  theme(strip.background = element_rect(colour = "black", fill="grey"),
        panel.border = element_rect(fill = NA),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 14)) 
print(p)

#------------------------------------------------------------------------------------------------
# Calculate change in %methylation status from WT for hotspots
#------------------------------------------------------------------------------------------------

mymat <- perc.meth.count.hotspots[,-c(NCOL(perc.meth.count.hotspots))]

for (k in 1:NROW(mymat)) {
  for (l in 2:NCOL(mymat)){
    if (mymat[k,l] < mymat[k,1]){
      #finding hotspots in WT that lose methylation (or become coldspots) in mutants, so subtracting %methylation status from WT
      mymat[k,l]=(mymat[k,1]-mymat[k,l])*100
    } else if (mymat[k,l] >= mymat[k,1]){
      mymat[k,l]=0
    }
  }
}
mymat

#------------------------------------------------------------------------------------------------------
#filtering for rows where there is at least 50% decrease in %methylation status in any mutant globally. 
#I used 10% cutoff for cmt3 and suvh456
#------------------------------------------------------------------------------------------------------
mymat.fltered <- mymat[apply(mymat>75, 1, any),]

#finding rows from original dataframe
orig.df <- perc.meth.count.hotspots[,-c(NCOL(perc.meth.count.hotspots))]
data_subset <- orig.df[which(rownames(orig.df) %in% rownames(mymat.fltered)),]
breaksList <- seq(from=min(data_subset), to=max(data_subset), length.out = 10)
#plotting original values
h <- pheatmap(data_subset, 
              cluster_cols=FALSE,
              #gaps_col = 1,
              breaks=breaksList,
              cluster_rows = FALSE,
              clustering_distance_rows = "euclidean", 
              #clustering_distance_cols = "euclidean", 
              clustering_method = "complete",
              gaps_col = c(1),
              show_rownames = FALSE,
              border_color=NA,
              main="%methylation status in WT and mutants",
              color = colorRampPalette(brewer.pal(n = 3, name = "YlOrRd"))(length(breaksList)))



