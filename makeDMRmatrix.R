
makeDMRmatrix <- function(context, chr, samplefile, input.dir, out.dir) {
  
  mylist <- list()
  
  # Read the sample file with filenames
  samplelist <- fread(samplefile, header=T)
  
  for (i in  1:length(context)){
    
    for (j in 1:length(chr)){
      selectlist <- list()
      cat(paste0("Extracting region files for ", chr[j], " ", context[i], "\n"), sep = "")
      
      # list all Region files in the input directory
      extractflist <- list.files(input.dir, pattern=paste0(chr[j],"_",context[i],".txt"), full.names=TRUE)
      
      if (length(extractflist) != 0){
        mynames <- gsub(paste0("_", chr[j], "_", context[i], ".txt$"), "", basename(extractflist))
        
        for (a1 in seq_along(mynames)){
          as <- samplelist[grepl(mynames[a1], samplelist$file),]
          as$full.path.regfile <- grep(paste0(mynames[a1], "_", chr[j], "_", context[i]), extractflist, value=TRUE)
          selectlist[[a1]] <- as
        }
        flist <- rbindlist(selectlist)
        print(flist)
        
        if (!is.null(flist$replicate)) {
          cat(paste0("Input data with replicates....creating unique sample names\n"), sep = "")
          flist$name <- paste0(flist$sample,"_", flist$replicate)
        } else {
          flist$name <- samplelist$sample 
        }
        
        cat(paste0("Now, constructing DMR matrix\n"), sep = "")
        # (column 6) state calls and (column 7) rc.meth.lvl for all samples; make 2 matrices
        merge.cols <- function(colm) {
          for (l in 1:length(colm)){
            extract <- lapply(flist$full.path.regfile, function(k){
              f <- fread(k, header=FALSE, skip=1, select=c(1, 2, 3, colm[l]))
              if (colm[l]==6) {
                f[,4] <- ifelse(f[,4] == "U", yes = 0, no = 1)
              }
              colnames(f)[4] <- basename(k)
              return(f)
            })
            df <- Reduce(function(x, y) {
              dplyr::inner_join(x, y, by=c("V1","V2","V3"))
            }, extract)
            
            # renaming file names with sample names
            for (a in 4:length(colnames(df))) {
              for (n in 1:length(flist$name)) {
                if (colnames(df)[a] == basename(flist$full.path.regfile)[n]) {
                  colnames(df)[a] = flist$name[n]
                }
              }
            }
            mylist[[l]] <- df
          }
          return(mylist)
        }
        
        #(column 6) state-calls and (column 7) rc.meth.lvl
        mydf <- merge.cols(colm=c(6, 7))
        
        # list containing the state calls
        status.collect <- mydf[[1]]
        # list containing the rcmethlvls
        rc.methlevel.collect <- mydf[[2]]
        
        names(status.collect)[1] <- "seqnames"
        names(status.collect)[2] <- "start"
        names(status.collect)[3] <- "end"
        
        names(rc.methlevel.collect)[1] <- "seqnames"
        names(rc.methlevel.collect)[2] <- "start"
        names(rc.methlevel.collect)[3] <- "end"
        
        cat("\n")
        cat(paste0("Writing output....\n"), sep = "")
        fwrite(x=status.collect, file=paste0(out.dir,"/", chr[j], "_", context[i],"_state-calls.txt"),
               quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
        fwrite(x=rc.methlevel.collect, file=paste0(out.dir,"/", chr[j], "_", context[i],"_rcmethlvl.txt"), 
               quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      } else{
        cat("Files do not exist\n")
      }
    }
  }
}
