
MethimputeRegTobedGraph <- function(regfile, out.dir) {
  cat(paste0("Reading Methimpute file: ", regfile, "\n"))
  fname <- fread(regfile, skip = 1, select = c("V1","V2","V3","V7"))
  name <- gsub(pattern = "\\.txt$", "", basename(regfile))
  fname$V2 <- format(fname$V2, scientific = FALSE) 
  fname$V3 <- format(fname$V3, scientific = FALSE)
  print(paste0("Writing to bedGraph format....",name))
  fwrite(x= fname, file=paste0(out.dir, "/", name,".bedGraph"), sep="\t", 
         quote=FALSE, row.names=FALSE, col.names=FALSE)
}
