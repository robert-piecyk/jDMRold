
modifiedExportMethylome <- function(model, out.dir, context, name) {
    data <- model$data
    final_dataset <- as(data, 'data.frame')
    final_dataset <- final_dataset[,c('seqnames','start','end','strand',
      'context','counts.methylated','counts.total',
      'posteriorMax','status','rc.meth.lvl')]

    # dropping columns
    drops <- c('width','strand','clusterlen','counts.methylated',
      'counts.total', 'distance', 'transitionContext', 'posteriorMeth','posteriorUnmeth')
    final_dataset <- final_dataset[ , !(names(final_dataset) %in% drops)]
    #------------------------------------------------------------------
    # convert full string into M/U/I
    final_dataset <- statusStringCheck(final_dataset)
    #-------------------------------------------------------------
    # take 4 digit of decimal value posteriorMax column
    final_dataset$posteriorMax <-floorDec(as.numeric(as.character(final_dataset$posteriorMax)),5)
    final_dataset$rc.meth.lvl <- floorDec(as.numeric(as.character(final_dataset$rc.meth.lvl)),5)
    final_dataset$seqnames <- as.character(final_dataset$seqnames)

    saveFile <- paste0(out.dir, name, "_", context, ".txt")
    fwrite(final_dataset, file = saveFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    return (final_dataset)
}

#--------------------------------------------------------------------------
makeRegionsImpute <- function(df, context, refRegion) {
    #regions file
    #data <- dget(regionFile)
    #data <- refRegion
    tmp_reg <- dget(refRegion)
    data <- as.data.frame(tmp_reg$reg.obs)
    #colnames(data)[which(names(data) == "cluster.size")] <- "cluster.length"
    
    #reference methimpute file G0
    ref_data <- fread(df, skip = 1, select = c("V1","V2","V3","V4","V5","V6"))
    ref_data <- ref_data %>% filter(ref_data$V1 != "M" & ref_data$V1 != "C")
    ref_data <- ref_data[which(ref_data$V4==context),]

    data_gr <- GRanges(seqnames=data$chr,
      ranges=IRanges(start=data$start, end=data$end),
      clusterlen=data$cluster.length,
      context=as.factor(context))

    ref_gr <- GRanges(seqnames=ref_data$V1,
      ranges=IRanges(start=ref_data$V2, width=1),
      context=as.factor(context),
      methylated=ref_data$V5,
      total=ref_data$V6)

    counts <- array(NA, dim=c(length(data_gr), 2),
      dimnames=list(NULL, c("methylated", "total")))

    overlaps <- findOverlaps(ref_gr, data_gr)


    overlaps.hits <- ref_gr[queryHits(overlaps)]

    methylated <- aggregate(overlaps.hits$methylated, list(subjectHits(overlaps)), FUN=sum)
    total <- aggregate(overlaps.hits$total, list(subjectHits(overlaps)), FUN=sum)
    #------------
    if (NROW(methylated) != NROW(counts) ){
    missingr <- which(!rownames(data.frame(data_gr)) %in% methylated$Group.1)

    for(item in seq_len(NROW((missingr)))){
        methylated <- rbind (c(missingr[item],0), methylated)
        methylated <- methylated[order(methylated$Group.1),]
        total <- rbind (c(missingr[item],0), total)
        total <- total[order(total$Group.1),]
    }
  }
  counts[,"methylated"] <- methylated$x
  counts[,"total"] <- total$x
  data_gr$counts <- counts
  return(data_gr)
}

#--------------------------------------------------------------------------
makeMethimpute<-function(df, context, refRegion, include.intermediate, probability, out.dir, name){
    methylome.data <- makeRegionsImpute(df, context, refRegion)
    quant.cutoff <- as.numeric(quantile(methylome.data$counts[,"total"], probs = c(0.96), na.rm=TRUE))
    distcor <- distanceCorrelation(methylome.data, distances=0:100)
    fit <- estimateTransDist(distcor)
    model <- callMethylation(data = methylome.data,
      transDist = fit$transDist,
      count.cutoff = quant.cutoff,
      max.iter = Inf,
      include.intermediate = include.intermediate,
      update = probability)
    methFile <- modifiedExportMethylome(model, out.dir, context, name)
    return(methFile)
}

#--------------------------------------------------------------------------
runMethimputeRegions <- function(Methfiles, Regionfiles, context, include.intermediate, probability, out.dir) {
  Methfiles <- list.files(Methfiles, pattern='\\.txt$', full.names = TRUE)
  for (i in 1:length(Methfiles)){
    for (j in 1:length(context)){
      methfn <- gsub(".*methylome_|\\.txt$", "", Methfiles[i])
      Regfiles <- list.files(Regionfiles, pattern=paste0("_",context[j],".Rdata"), full.names = TRUE)
      for (k in 1:length(Regfiles)){

        tmp <- gsub(".*Arabidopsis_regions_|\\.Rdata$", "", Regfiles[k])
        chr <- gsub(paste0("_", context[j]), "", tmp)

        name <- paste0(methfn, "_", chr)
        cat(paste0("Running for ", methfn, " ", context[j], " ", chr,"\n"), sep = "")
        makeMethimpute(
          df=Methfiles[i],
          context=context[j],
          refRegion=Regfiles[k],
          include.intermediate=include.intermediate, 
          probability=probability,
          out.dir=out.dir,
          name=name)
      }
    }
  }
}
