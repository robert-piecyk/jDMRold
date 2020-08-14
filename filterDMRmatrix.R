
filterDMRmatrix <- function(epiMAF.cutoff, replicate.consensus, data.dir) {
  
  list.status <- list.files(data.dir, pattern="_state-calls.txt", full.names=TRUE)
  
  for (i in seq_along(list.status)){
    
    fname <- gsub("_state-calls.txt", "", basename(list.status[i]))
    
    cat(paste0("Running DMR matrix for ", fname, "\n"), sep = "")
    cat("\n")
    #----------------------------------------------
    # Removing non-polymorphic/unchanged patterns
    #----------------------------------------------
    
    status.collect <- fread(list.status[i], header=T)
    rc.methlevel.collect <- fread(paste0(data.dir, "/", fname, "_rcmethlvl.txt"), header=T)
    
    cat(paste0("Removing non-polymorphic patterns......\n"))
    
    index <- which(rowSums(status.collect[,4:NCOL(status.collect)]) != 0 & 
                     rowSums(status.collect[,4:NCOL(status.collect)]) != NCOL(status.collect)-3)
    
    status.collect <- status.collect[index,]
    rc.methlevel.collect <- rc.methlevel.collect[index,]
    
    #----------------------------------------------
    # Filtering out regions with epiMAF < Minor Allele Frequency
    #----------------------------------------------
    
    if (!is.null(epiMAF.cutoff)) {
      
      cat(paste0("Filtering for epiMAF......\n"))
      
      status.collect$epiMAF <- 0
      
      for (i1 in 1:NROW(status.collect)){
        mypattern <- unlist(status.collect[i1, 4:(NCOL(status.collect)-1)])
        mycount <- table(mypattern)
        epiMAF <- min(mycount)/length(mypattern)
        status.collect$epiMAF[i1] <- floorDec(as.numeric(as.character(epiMAF)),5)
      }
      status.collect <- status.collect[which(status.collect$epiMAF < epiMAF.cutoff),]
    }
    
    #----------------------------------------------
    # Retaining samples based on replicate.consensus
    #----------------------------------------------
    
    if (!is.null(replicate.consensus)) {
      
      if (!is.null(status.collect$epiMAF)){
        status.collect <- status.collect[,-c("epiMAF")]
      }
      
      # deducing replicates info
      mycol <- names(status.collect)[4:NCOL(status.collect)]
      sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))
      colnames(sampleinfo) <- c("sample", "replicate")
      
      cat(paste0("Filtering for replicate consensus.....\n"))
      cat("\n")
      
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
        }
        #print(df.bind)
        if (sum(df.bind$count)==0)
          dt <- rbind(dt, status.collect[x,])
      })
      status.collect <- q[-which(sapply(q, is.null))]
      status.collect <- rbindlist(status.collect)
    }
    
    df.rc.methlevel.collect <- rc.methlevel.collect %>% semi_join(status.collect, by=c("seqnames","start","end"))
    
    fwrite(x=status.collect, file=paste0(data.dir, "/", fname,"_state-calls-filtered.txt"), quote=FALSE, 
           row.names=FALSE, col.names=TRUE, sep="\t")
    fwrite(x=df.rc.methlevel.collect, file=paste0(data.dir,"/", fname,"_rcmethlvl-filtered.txt"), quote=FALSE, 
           row.names=FALSE, col.names=TRUE, sep="\t")
  }
}
