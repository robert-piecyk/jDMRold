
filterDMRmatrix <- function(context, chr, epiMAF.cutoff, replicate.consensus, out.dir) {
  
  for (i in  1:length(context)){
    for (j in 1:length(chr)){
      cat(paste0("Running for ", chr[j], "_", context[i], ":\n"), sep = "")
      cat("\n")
      
      status.collect <- fread(paste0(out.dir, chr[j], "_", context[i], "_state-calls.txt", sep=""))
      rc.methlevel.collect <- fread(paste0(out.dir, chr[j], "_", context[i], "_rcmethlvl.txt", sep=""))
      
      #----------------------------------------------
      # For any kind of data. Removing non-polymorphic/unchanged patterns
      
      cat(paste0("Removing non-polymorphic patterns......\n"))
      cat("\n")
      
      index <- which(rowSums(status.collect[,4:NCOL(status.collect)]) != 0 & 
                       rowSums(status.collect[,4:NCOL(status.collect)]) != NCOL(status.collect)-3)
      
      status.collect <- status.collect[index,]
      rc.methlevel.collect <- rc.methlevel.collect[index,]
      
      #----------------------------------------------
      # Filtering out regions for Minor Allele Frequency (lower than epiMAF cutoff)
      
      if (!is.null(epiMAF.cutoff)) {
        
        cat(paste0("Filtering for epiMAF......\n"))
        cat("\n")
        
        status.collect$epiMAF <- 0
        
        for (i1 in 1:NROW(status.collect)){
          mypattern <- unlist(status.collect[i1, 4:(NCOL(status.collect)-1)])
          mycount <- table(mypattern)
          epiMAF <- min(mycount)/length(mypattern)
          status.collect$epiMAF[i1] <- floorDec(as.numeric(as.character(epiMAF)),5)
        }
        status.collect <- status.collect[which(status.collect$epiMAF < epiMAF.cutoff),]
      } 
      
      #-----------------------------------------------
      # Retaining samples based on replicate.consensus
      
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
            pattern.vals  <- df.bind$mypattern[df.bind$sample==m]
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
      
      fwrite(x=status.collect, file=paste0(out.dir, chr[j], "_", context[i],"_state-calls-filtered.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=df.rc.methlevel.collect, file=paste0(out.dir, chr[j], "_", context[i],"_rcmethlvl-filtered.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
    }
  }
}


