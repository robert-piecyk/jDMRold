
filterReplicateConsensus <- function(status.collect, rc.methlevel.collect, replicate.consensus, gap=1){
  
  pb1 <- txtProgressBar(min = 1, max = NROW(status.collect), char = "=", style = 3, file = "")
  
  if (!is.null(status.collect$epiMAF)){
    status.collect <- status.collect[,-c("epiMAF")]
  }
  # deducing replicates info
  mycol <- names(status.collect)[4:NCOL(status.collect)]
  sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))
  colnames(sampleinfo) <- c("sample", "replicate")
  
  dt <- data.frame()
  q <- lapply(1:NROW(status.collect), function(x){
    mypattern <- unlist(status.collect[x, 4:NCOL(status.collect)])
    df.bind <- cbind(sampleinfo, mypattern)
    for (m in unique(df.bind$sample)){
      rval <- replicate.consensus * length(df.bind$mypattern[df.bind$sample==m])
      pattern.vals <- df.bind$mypattern[df.bind$sample==m]
      tt <- table(pattern.vals)
      if (max(tt) >= rval)
        df.bind$count[df.bind$sample==m] <- 0
      else
        df.bind$count[df.bind$sample==m] <- 1
      
      Sys.sleep(1/NROW(status.collect))
      setTxtProgressBar(pb1, x)
    }
    close(pb1)
    #print(df.bind)
    if (sum(df.bind$count)==0)
      dt <- rbind(dt, status.collect[x,])
  })
  status.collect <- q[-which(sapply(q, is.null))]
  df.status.collect <- rbindlist(status.collect)
  if (NROW(df.status.collect) !=0){
  df.rc.methlevel.collect <- rc.methlevel.collect %>% semi_join(df.status.collect, by=c("seqnames","start","end"))
  return(list(df.status.collect, df.rc.methlevel.collect))
  } else {
    message("\nEmpty dataframe. Nothing to write!")
    return(NULL)
  }
}


filterEpiMAF <- function(mat1, mat2, epiMAF){
  
  pb2 <- txtProgressBar(min = 1, max = NROW(mat1), char = "=", style = 3, file = "")
  mat1$epiMAF <- 0
  
  for (i1 in 1:NROW(mat1)){
    mypattern <- unlist(mat1[i1, 4:(NCOL(mat1)-1)])
    mycount <- table(mypattern)
    epiMAF.out <- min(mycount)/length(mypattern)
    mat1$epiMAF[i1] <- floorDec(as.numeric(as.character(epiMAF.out)),5)
    
    Sys.sleep(1/NROW(mat1))
    setTxtProgressBar(pb2, i1)
  }
  close(pb2)

  df.status.collect <- mat1[which(mat1$epiMAF < epiMAF),]
  if (NROW(df.status.collect) !=0){
    df.rc.methlevel.collect <- mat2 %>% semi_join(df.status.collect, by=c("seqnames","start","end"))
    return(list(df.status.collect, df.rc.methlevel.collect))
  } else {
    message("\nEmpty dataframe. Nothing to write!")
    return(NULL)
  }
}

merge.bins <- function(rcmethlvl, statecalls, gap){
  mylist <- list()
  result <- list()
  
  matrix1 <- as.data.frame(statecalls)
  matrix2 <- as.data.frame(rcmethlvl)
  
  # extract unique pattern
  extract.pattern <- unique(matrix1[,(4:NCOL(matrix1))])
  extract.pattern$pattern <- apply(extract.pattern[,c(1:NCOL(extract.pattern))], 1, paste, collapse="")
  
  # the state call matrix
  gr1 <- GRanges(seqnames=matrix1$seqnames, ranges=IRanges(start=matrix1$start, end=matrix1$end))
  values(gr1) <- cbind(values(gr1), pattern=apply(matrix1[,c(4:NCOL(matrix1))], 1, paste, collapse=""))
  
  # the rcmethlvl matrix, also add the state-call pattern
  gr2 <- GRanges(seqnames=matrix2$seqnames, ranges=IRanges(start=matrix2$start, end=matrix2$end))
  values(gr2) <- cbind(values(gr2), values(gr1), DataFrame(matrix2[,c(4:NCOL(matrix2))]))
  
  # this is for the state-calls: collapse overlapping bins if pattern is same
  grl_reduce <- unlist(GenomicRanges::reduce(split(gr1, gr1$pattern)))
  result <- sort(grl_reduce)
  result$pattern <- names(result)
  result <- data.frame(result)
  final.status.collect <- result %>% left_join(extract.pattern, by=c("pattern"))
  
  # this is for the rcmethlvl: collapse bins and take average of the bins
  mycols <- colnames(matrix1)[4:NCOL(matrix1)]
  fn = function(u){
    out = GenomicRanges::reduce(u)
    for (x in 1:length(mycols)){
      eval(parse(text=paste0("out$", mycols[x], " = mean(u$", mycols[x], ")")))
    }
    return(out)
  }
  
  message("\nNow, Merging overlapping and consecutive bins...\n")
  # split rcmthlvl matrix based on pattterns
  grl <- split(gr2, gr2$pattern)
  
  pb3 <- txtProgressBar(min = 1, max = length(grl), char = "=", style = 3, file = "")
  
  for (x in 1:length(grl)){
    a <- data.frame(grl[[x]])
    a1 <- a %>% arrange(pattern, start) %>% group_by(pattern) %>% 
      mutate(indx = cumsum(start > lag(end, default = start[1]) + gap))
    a1.gr <- makeGRangesFromDataFrame(a1, keep.extra.columns=TRUE)
    df <- lapply(split(a1.gr, a1.gr$indx), fn)
    mylist[[x]] <- df
    
    Sys.sleep(1/length(grl))
    setTxtProgressBar(pb3, x)
  }
  close(pb3)
  
  f.df <- unlist(mylist)
  f.df <- do.call(rbind, lapply(f.df, data.frame))
  f.df <- f.df[order(f.df[,1], f.df[,2]),]
  f.df[,(6:NCOL(f.df))] <- lapply(f.df[,(6:NCOL(f.df))], function(xy){ floorDec(xy,5)})
  
  final.status.collect <- subset(final.status.collect, select = -c(strand, pattern))
  final.rcmethlvl.collect <- subset(f.df, select = -c(strand))
  return(list(final.status.collect, final.rcmethlvl.collect))
}

export.out <- function(out.rcmethlvl, out.statecalls, context, out.name1, out.name2, data.out){
  fwrite(x=out.statecalls, file=paste0(data.out, "/", context, "_", out.name1, ".txt"), 
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  fwrite(x=out.rcmethlvl, file=paste0(data.out, "/", context, "_", out.name2, ".txt"), 
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

filterDMRmatrix <- function(epiMAF.cutoff, replicate.consensus, gridDMR, data.dir) {
  
  list.status <- list.files(data.dir, pattern="_StateCalls.txt", full.names=TRUE)
  if (length(list.status) != 0){
    for (i in seq_along(list.status)){
      context <- gsub("_StateCalls.txt", "", basename(list.status[i]))
      cat("\n")
      cat(paste0("Running DMR matrix for ", context, "\n"), sep = "")
      cat("\n")
      
      #----------------------------------------------
      # Removing non-polymorphic/unchanged patterns
      #----------------------------------------------
      status.collect <- fread(list.status[i], header=T)
      rc.methlevel.collect <- fread(paste0(data.dir, "/", context, "_rcMethlvl.txt"), header=T)
      
      cat(paste0("Removing non-polymorphic patterns...\n"))
      
      index <- which(rowSums(status.collect[,4:NCOL(status.collect)]) != 0 & 
                       rowSums(status.collect[,4:NCOL(status.collect)]) != NCOL(status.collect)-3)
      
      status.collect <- status.collect[index,]
      rc.methlevel.collect <- rc.methlevel.collect[index,]
      
      if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
        message("Both, epiMAF and replicate consensus set to NULL")
        out1=status.collect
        out2=rc.methlevel.collect
        export.out(out.statecalls=out1,
                   out.rcmethlvl=out2,
                   context=context,
                   out.name1="StateCalls-filtered",
                   out.name2="rcMethlvl-filtered",
                   data.out=data.dir)
      }
      #----------------------------------------------
      # Optional. Filtering out regions with epiMAF < Minor Allele Frequency
      #----------------------------------------------
      if (!is.null(epiMAF.cutoff)) {
        cat(paste0("Filtering for epiMAF: ", epiMAF.cutoff, "\n"))
        cat("\n")
        mydf <- filterEpiMAF(mat1=status.collect, mat2=rc.methlevel.collect, epiMAF=epiMAF.cutoff)
        if (!is.null(mydf)){
          # For Population data remove the epiMAF column
          out1=mydf[[1]][,-c("epiMAF")]
          out2=mydf[[2]]
          export.out(out.statecalls=mydf[[1]],
                     out.rcmethlvl=mydf[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
      } 
      
      #----------------------------------------------
      # Optional. Retaining samples based on replicate.consensus
      #----------------------------------------------
      if (!is.null(replicate.consensus)) {
        cat(paste0("Filtering for replicate consensus...\n"))
        cat("\n")
        mydf <- filterReplicateConsensus(status.collect, rc.methlevel.collect, replicate.consensus)
        if (!is.null(mydf)){
          out1=mydf[[1]]
          out2=mydf[[2]]
          export.out(out.statecalls=out1,
                     out.rcmethlvl=out2,
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
      }
      
      #----------------------------------------------
      # Merging bins
      #----------------------------------------------
      if (gridDMR==TRUE) {
        if (exists("out1") && exists("out2")){
          out <- merge.bins(statecalls=out1, rcmethlvl=out2, gap=1)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered-merged",
                     out.name2="rcMethlvl-filtered-merged",
                     data.out=data.dir)
        } 
      } else {
        message("\ngrid DMR set to FALSE")
        }
    }
  } else {
    message("\nDMR matrix files do not exist!")
  }
}