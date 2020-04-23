
# set number of cores based on cpu, number of jobs
numCore <- function(NumFiles){
  no_cores <- detectCores()
  if (no_cores > NumFiles){
    no_cores <- NumFiles
  }else{
    no_cores <- no_cores - 1
  }
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  return(cl)
}

# take 5 digit of decimal value and cut the numbers
floorDec <- function(valParm ,x){
  y <- function(x, level=1) round(x - 5*10^(-level-1), level)
  res <-y(as.numeric(valParm),x)
  return(res)
}


# Name of data files
getNames <- function(fileName){
  tmpName <- gsub(pattern = "\\methylome_", "", basename(fileName))
  fileName <- gsub(pattern = "\\.txt", "",tmpName)
  return(fileName)
}


# Convert all String into U/I/M if there are full text.
statusStringCheck <-  function(file_A){
  list_status <- c("Unmethylated", "Intermediate", "Methylated")
  strTocheckFileA <- utils::head(file_A$status[1])
  if (strTocheckFileA %in% list_status) {
    file_A$status <- str_replace_all(file_A$status,
      pattern = "Unmethylated", replacement = "U")
    file_A$status <- str_replace_all(file_A$status,
      pattern = "Intermediate", replacement = "I")
    file_A$status <- str_replace_all(file_A$status,
      pattern = "Methylated", replacement = "M")
  }
  return(file_A)
}


# saving results- DMR pattern
saveResult<-function(finalDF, name, cytosine, out.dir, current){
  saveFile <- paste0(out.dir, name, cytosine,"_", current, ".csv")
  fwrite(finalDF, file = saveFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  message("DMR Results saved in: ", saveFile, "\n")

}
# saving results- same pattern
saveSameCyt<-function(samePattern, name, cytosine, out.dir, current){
  saveSamepattern <- paste0(out.dir, name, cytosine, "_", current, ".csv")
  fwrite(samePattern,file = saveSamepattern ,quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
  message("Same Cytosines: ", saveSamepattern, "\n")

}
