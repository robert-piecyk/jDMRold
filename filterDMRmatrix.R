
filterDMRmatrix <- function(context, chr, out.dir) {

  for (i in  1:length(context)){
    for (j in 1:length(chr)){
      cat(paste0("Running for ", chr[j], "_", context[i], "\n"), sep = "")
      
      status.collect <- fread(paste0(out.dir, chr[j], "_", context[i], "_state-calls.txt", sep=""))
      rc.methlevel.collect <- fread(paste0(out.dir, chr[j], "_", context[i], "_rcmethlvl.txt", sep=""))
     
      # deducing replicates info
      mycol <- names(status.collect)[4:NCOL(status.collect)]
      sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))
      colnames(sampleinfo) <- c("sample", "replicate")
      
      # Retaining replicate consensus among samples
      dt <- data.frame()
      q <- lapply(1:NROW(status.collect), function(x){
        mypattern <- unlist(status.collect[x, 4:NCOL(status.collect)])
        df.bind <- cbind(sampleinfo, mypattern)
      
        for (m in unique(df.bind$sample)){
          if (length(unique(df.bind$mypattern[df.bind$sample==m]))==1)
            df.bind$count[df.bind$sample==m] <- 0
          else
            df.bind$count[df.bind$sample==m] <- 1
        }
        #print(df.bind)
        if (sum(df.bind$count)==0)
          dt <- rbind(dt, status.collect[x,])
      })
      df.status.collect <- q[-which(sapply(q, is.null))]
      df.status.collect <- rbindlist(df.status.collect)
      
      df.rc.methlevel.collect <- rc.methlevel.collect %>% semi_join(df.status.collect, by=c("seqnames","start","end"))
      
      fwrite(x=df.status.collect, file=paste0(out.dir, chr[j], "_", context[i],"_state-calls-filtered.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=df.rc.methlevel.collect, file=paste0(out.dir, chr[j], "_", context[i],"_rcmethlvl-filtered.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
    }
  }
}


