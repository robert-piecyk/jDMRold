
makeDMRmatrix <- function(context, chr, samplefile, input.dir, out.dir) {
  
  mylist <- list()
  
  # Read the sample file with filnames
  samplelist <- fread(samplefile, header=T)
  
  for (i in  1:length(context)){
    for (j in 1:length(chr)){
      cat(paste0("Extracting region files for ", chr[j], " ", context[i], "\n"), sep = "")
      extractflist <- list.files(input.dir, 
                             pattern=paste0(chr[j],"_",context[i],".txt"), 
                             full.names=TRUE)
      
      # Matching filenames
      samplelist$full.path <- NA
      for (a in 1:length(extractflist)) {
        for (n in 1:length(samplelist$file)) {
          pat1 <- gsub(".*methylome_|\\_All.txt$", "", samplelist$file[n])
          pat2 <- gsub(paste0("_",chr[j],"_",context[i],".txt"), "", basename(extractflist[a]))
          if (pat1 == pat2) {
            samplelist$full.path[n] <- extractflist[a]
          }
        }
      }
  
      cat("\n")
      #extractp <- paste0(chr[j],"_",context[i])
      #samplelist <- samplelist[grep(input.dir, samplelist$file[i]),]
      
      if (!is.null(samplelist$replicate)) {
        cat(paste0("Input data with replicates....creating unique sample names\n"), sep = "")
        samplelist$name <- paste0(samplelist$sample,"_", samplelist$replicate)
      } else {
        samplelist$name <- samplelist$sample }
      
      print (samplelist[,c("file","full.path","name")])
      
      cat(paste0("Now, constructing DMR matrix\n"), sep = "")
      # (column 6) state calls and (column 7) rc.meth.lvl for all samples and make 2 matrices
      merge.cols <- function(colm) {
        for (l in 1:length(colm)){
          extract <- lapply(samplelist$full.path, function(k){
            f <- fread(k, header=FALSE, skip=1, 
                       select=c(1, 2, 3, colm[l]))
            if (colm[l]==6) {
              f[,4] <- ifelse(f[,4] == "U", yes = 0, no = 1)
            }
            colnames(f)[4] <- basename(k)
            return(f)
          })
          df <- Reduce(function(x, y) {
            dplyr::inner_join(x, y, by=c("V1","V2","V3"))
          }, extract)
          
          # renaming file name with a sample name
          for (a in 4:length(colnames(df))) {
            for (n in 1:length(samplelist$name)) {
              if (colnames(df)[a] == basename(samplelist$full.path)[n]) {
                colnames(df)[a] = samplelist$name[n]
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
    }
  }
}
