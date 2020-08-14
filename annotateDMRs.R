
#-----------------------------------------------------------------------------------
# output DMR overlaps in gff3 format
#-----------------------------------------------------------------------------------
gff3.out <- function(annotation, grangesObj, gff, name, out.dir) {
  getgff3 <- lapply(annotation, function(x){
    idx <- which(elementMetadata(gff)[,"type"] == x)
    gff <- gff[idx,]
    hits <- findOverlaps(grangesObj, gff, ignore.strand=FALSE)
    gr.matched <-grangesObj[queryHits(hits)]
    mcols(gr.matched) <- cbind.data.frame(mcols(gr.matched), mcols(gff[subjectHits(hits)]))
    values(gr.matched) <- cbind(values(gr.matched), region="DMR")  
    names(elementMetadata(gr.matched))[names(elementMetadata(gr.matched)) == "type"] <- "annotation"
    return(gr.matched)
  })
  export.gff(do.call(c, getgff3), paste0(out.dir,"/", name, "_annotation.gff3"), version="3")
}

#-----------------------------------------------------------------------------------
# DMR counts
#-----------------------------------------------------------------------------------
annotateDMRs <- function(annotation, gff, gff3.out, file.list, out.dir) {
  
  anno.list <- list()
  final.df <- list()
  
  for (i in 1:length(file.list)){
    
    cat(paste0("Running file ", file.list[i], "\n"), sep = "")
    
    # standard input files are *state-calls-filtered.txt
    if (grepl("\\_state-calls-filtered.txt$", basename(file.list[i]))) {
      tmp.name <- gsub("\\_state-calls-filtered.txt$", "", basename(file.list[i]))
    } else {
      tmp.name <- gsub("\\.txt$", "", basename(file.list[i]))
    }
    
    file <- fread(file.list[i], skip=1, select=c(1,2,3))
    
    #if (is.null(file[,c("seqnames","start","end")])) {
    #  stop("Please supply input files with 3 columns: seqnames, start, end")
    #}
    
    gr <- GRanges(seqnames=file$V1, ranges=IRanges(start=file$V2, end=file$V3))
    
    if (gff3.out==TRUE) {
      gff3.out(annotation=annotation, gff=gff, grangesObj=gr, name=tmp.name, out.dir=out.dir)
    }
    
    # count the DMR overlaps; the output can be used to make a barplot or pie-chart
    getAnno <- lapply(annotation, function(x){
      idx <- which(elementMetadata(gff)[,"type"] == x)
      gff <- gff[idx,]
      hits <- findOverlaps(gr, gff, ignore.strand=FALSE)
      myranges <- subsetByOverlaps(gr, gff)
      
      mcols(myranges)$id <- CharacterList(split(gff$ID[subjectHits(hits)], queryHits(hits)))
      mcols(myranges)$type <- CharacterList(split(gff$type[subjectHits(hits)], queryHits(hits)))
      mcols(myranges)$chr <- CharacterList(split(seqnames(gff)[subjectHits(hits)], queryHits(hits)))
      mcols(myranges)$coord <- CharacterList(split(ranges(gff)[subjectHits(hits)], queryHits(hits)))
      mcols(myranges)$anno.str <- CharacterList(split(strand(gff)[subjectHits(hits)], queryHits(hits)))
      df <- as(myranges, "data.frame")
      df <- df %>% mutate(id = strsplit(as.character(id), ","),
                          type = strsplit(as.character(type), ","), 
                          chr = strsplit(as.character(chr), ","),
                          coord = strsplit(as.character(coord), ","),
                          anno.str = strsplit(as.character(anno.str), ",")) %>% 
        tidyr::unnest(c(id, type, chr, coord, anno.str))
      cleandf <- data.frame(lapply(df, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
      cleandf$Annotation.coord <- apply(cleandf[,c("chr","coord","anno.str")], 1, paste, collapse=":")
      cleandf <- subset(cleandf, select = -c(chr, coord, anno.str))
      return(cleandf)
    })
    
    d <- rbindlist(getAnno)
    out <- d %>%
      dplyr::group_by(seqnames, start, end) %>%
      dplyr::summarize(
        type = paste(type, collapse=","),
        id = paste(id, collapse=","),
        Annotation.coord = paste(Annotation.coord, collapse=","), .groups = 'drop'
      )
    # unique annotation overlaps
    out.1 <- out[which(sapply(strsplit(out$type,','), uniqueN)==1),]
    
    # multiple annotation overlaps
    mx <- out[which(sapply(strsplit(out$type,','), uniqueN)>=2),]
    
    if (nrow(d) != 0){
      for (k in seq_along(annotation)){
        anno.list[[k]] <- NROW(out.1[grep(annotation[k], out.1$type),])
        names(anno.list)[[k]] <- annotation[k]
      }
      df.1 <- do.call(cbind, anno.list)
      #print(anno.list)
      
      df <- cbind(sample=tmp.name, total.DMRs=NROW(gr), df.1, multiple.overlaps=NROW(mx))
      final.df <- rbind(final.df, data.frame(df))
      
      fwrite(x=out, file=paste0(out.dir, "/", tmp.name, "_annotation.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
    }
    fwrite(x=final.df, file=paste0(out.dir, "/DMR-counts.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  }
}

