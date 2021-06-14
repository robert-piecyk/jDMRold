
runMethimputeRegions <- function(samplefiles, Regionfiles, genome, context, out.dir, include.intermediate=FALSE, mincov=0, nCytosines=0) {
  df.obs <- list()
  df.sim <- list()
  merge.list <- vector(mode="list")
  
  # Read the sample file with filenames and file paths
  filelist <- fread(samplefiles, header=TRUE)

  for (j in 1:length(context)){
    Regfiles <- list.files(Regionfiles, pattern=paste0("_", context[j], ".Rdata"), full.names = TRUE)
    #print(Regfiles)
    cat(paste0("Reading Region files. Merging individual chr data for context ", context[j], " ...\n"), sep="")
    cat("\n")
    for (k1 in 1:length(Regfiles)){
      f.file <- dget(Regfiles[k1])
      if (NROW(f.file$reg.obs)==0) {
        cat(paste0("Empty file ", basename(Regfiles[k1]), " ...\n"), sep="")
      } else {
        df.obs[[k1]] <- as.data.frame(f.file$reg.obs)
        df.sim[[k1]] <- as.data.frame(f.file$reg.sim)
      }
    }
    outlist <- list(reg.obs=do.call(rbind,df.obs),
                    reg.sim=do.call(rbind,df.sim),
                    context=context[j])
    
    regMerged <- paste(out.dir, "/",genome,"_Merged_Regions_", context[j], ".Rdata", sep="")
    dput(outlist, regMerged)
    
    for (k2 in 1:length(filelist$file)){
      methfn <- gsub(".*methylome_|\\_All.txt$", "", filelist$file[k2])
      cat(paste0("Now running file: ", methfn, " for context ", context[j], " ...\n"), sep="")
      regions.out <- makeMethimpute(
        df=filelist$file[k2],
        context=context[j],
        refRegion=regMerged,
        fit.plot=FALSE,
        include.intermediate=include.intermediate, 
        probability="constrained",
        out.dir=out.dir,
        fit.name=paste0(methfn, "_", context[j]),
        name=methfn,
        nCytosines=nCytosines,
        mincov=mincov)
    }
  }
}

binGenome <- function(fasta, win, step, genome, out.dir){
  val <- c()
  cat("Extracting chromosomes\n")
  chr.ext <- list.files(fasta, pattern="fa|fasta.gz$|fa.gz$", include.dirs = FALSE, full.names=TRUE)
  if (length(chr.ext)==0) {
    stop ("Empty folder!")
  } else {
    for (x in seq_along(chr.ext)){
      f.name <- gsub("Arabidopsis_thaliana.TAIR10.dna.chromosome.|\\.fa|\\.fa.gz|\\.fasta.gz$", "", basename(chr.ext[x]))
      f <- seqinr::read.fasta(chr.ext[x])
      cat(paste0("Chr: ", f.name, " \n"))
      if (length(f)>1) {
        stop ("Exiting... Multi FASTA file detected! Please, supply individual fasta files\n") }
      else {
        val[x] <- getLength(f)
        names(val)[x] <- f.name
      } 
    }
  }
  
  #gr <- GRanges(seqnames=gsub("chr", "", names(chrfile)), ranges=IRanges(start=1, end=chrfile))
  gr <- GRanges(seqnames=names(val), ranges=IRanges(start=1, end=val))
  
  binned.g <- slidingWindows(gr, width = win, step = step)
  d <- data.frame(unlist(binned.g))
  names(d)[1] <- "chr"
  names(d)[2] <- "start"
  names(d)[3] <- "end"
  names(d)[4] <- "cluster.length"
  cat(paste0("Binning genome with windows of: ", win, " bp and step-size of: ", step, " bp\n"), sep = "")
  new <- list(d)
  names(new) <- "reg.obs"
  out.name <- paste0(out.dir, "/", genome,"_Win", win, "_Step", step, ".Rdata", sep="")
  dput(new, out.name)
}

runMethimputeGrid <- function(out.dir, fasta, win, step, genome, samplefiles, context, mincov, include.intermediate=FALSE, nCytosines){

  binGenome(fasta, win, step,genome, out.dir)
  merge.list <- vector(mode="list") 
  filelist <- fread(samplefiles, header=TRUE)
  for (j in 1:length(context)){
    for (i in 1:length(filelist$file)){
      methfn <- gsub(".*methylome_|\\.txt|_All.txt$", "", filelist$file[i])
      cat(paste0("Running file: ",methfn," for context: ",context[j],"\n"), sep = "")
      grid.out <- makeMethimpute(
        df=filelist$file[i],
        context=context[j],
        refRegion=paste0(out.dir,"/",genome,"_Win",win,"_Step",step,".Rdata",sep=""),
        fit.plot=FALSE,
        include.intermediate=include.intermediate, 
        probability="constrained",
        out.dir=out.dir,
        fit.name=paste0(methfn, "_", context[j]),
        name=methfn,
        nCytosines=nCytosines,
        mincov=mincov)
    }
  }
}
