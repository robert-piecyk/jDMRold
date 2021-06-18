#'
#' @param gff
#' @importFrom  rtracklayer import.gff3
#' @export
#' @return merge all supplied gff3 annotations into one
#merge the input gff3 files into one
gff3.in <- function(gff){
  input.gff <- lapply(gff, function(x){
    import.gff3(x, colnames=c("type", "ID"))
  })
  merged.gff <- do.call(c, input.gff)
  return(merged.gff)
}

# Check Annotation levels here. Supply annotation terms for e.g genes, TEs
#levels(elementMetadata(merged.gff)[,"type"])
#available annotations
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"

#' @param gff
#' @param annotation
#' @param grangesObj
#' @param name
#' @param out.dir output directory
#' @import GenomicRanges
#' @importFrom  rtracklayer export.gff
#' @export
#' @return export output files in gff3 format
#output annotated gff3 files
gff3.out <- function(annotation, grangesObj, gff, name, out.dir) {
  getgff3 <- lapply(annotation, function(x){
    idx <- which(elementMetadata(gff)[,"type"] == x)
    gff <- gff[idx,]
    hits <- findOverlaps(grangesObj, gff, ignore.strand=FALSE)
    gr.matched <- grangesObj[queryHits(hits)]
    if (NROW(gr.matched) != 0){
      mcols(gr.matched) <- cbind.data.frame(mcols(gr.matched), mcols(gff[subjectHits(hits)]))
      values(gr.matched) <- cbind(values(gr.matched), region="DMR")
      names(elementMetadata(gr.matched))[names(elementMetadata(gr.matched)) == "type"] <- "annotation"
    }
    return(gr.matched)
  })
  export.gff(do.call(c, getgff3), paste0(out.dir,"/", name, "_annotation.gff3"), version="3")
}

#' @inheritParams gff3.in
#' @param gff
#' @param annotation
#' @param grangesObj
#' @param name
#' @param out.dir output directory
#' @param getAnno
#' @param mygff
#' @param mygr
#' @param annotation
#' @param gff.files
#' @param gff3.out
#' @param input.dir
#' @import  GenomicRanges
#' @export
#' @return annotated list
#extract annotated regions
annotate <- function(getAnno, mygff, mygr){
  lapply(getAnno, function(x){
    idx <- which(elementMetadata(mygff)[,"type"] == x)
    mygff <- mygff[idx,]
    hits <- findOverlaps(mygr, mygff, ignore.strand=FALSE)
    myranges <- subsetByOverlaps(mygr, mygff)

    mcols(myranges)$id <- CharacterList(split(mygff$ID[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$type <- CharacterList(split(mygff$type[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$chr <- CharacterList(split(seqnames(mygff)[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$coord <- CharacterList(split(ranges(mygff)[subjectHits(hits)], queryHits(hits)))
    mcols(myranges)$anno.str <- CharacterList(split(strand(mygff)[subjectHits(hits)], queryHits(hits)))
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
}

#' Annotate DMRs
#'
#' This function takes gff3 files as input and outputs annotated DMRs in text and gff3 format. Additionally, a DMR count table is generated.
#'
#' @inheritParams gff3.in
#' @inheritParams gff3.out
#' @inheritParams annotate
#' @param out.dir output directory
#' @param annotation annotation terms used to annotate DMRs
#' @param gff.files multiple gff3 files can be supplied as a vector
#' @param gff3.out a logical specifying whether output annotated files in gff3 format
#' @param input.dir input directory containing filtered DMR matrix/matrices. Ideally any file containing 3 columns i.e (chr, start, stop) can be supplied.
#' @import GenomicRanges
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  dplyr group_by
#' @importFrom  dplyr summarize
#' @export
#' @return output files containing annotated DMRs and DMR counts table.

annotateDMRs <- function(annotation, gff.files, gff3.out, input.dir, out.dir) {
  anno.list <- list(); final.df <- list()
  file.list <- list.files(input.dir, pattern="*.txt", full.names = TRUE)
  for (i in 1:length(file.list)){

    cat(paste0("Running file ", file.list[i], "\n"), sep = "")
    tmp.name <- gsub("\\.txt$", "", basename(file.list[i]))
    file <- fread(file.list[i], skip=1, select=c(1,2,3))
    gr <- GRanges(seqnames=file$V1, ranges=IRanges(start=file$V2, end=file$V3))
    gff=gff3.in(gff.files)

    if (gff3.out==TRUE) {
      gff3.out(annotation=annotation,
               gff=gff,
               grangesObj=gr,
               name=tmp.name,
               out.dir=out.dir)
    }

    d <- rbindlist(annotate(getAnno=annotation, mygff=gff, mygr=gr))
    if (nrow(d) != 0){
      out <- d %>%
        dplyr::group_by(seqnames, start, end) %>%
        dplyr::summarize(
          type = paste(type, collapse=","),
          id = paste(id, collapse=","),
          Annotation.coord = paste(Annotation.coord, collapse=",")
          #, .groups = 'drop'
          ) %>% as.data.frame()

      if (NROW(out)!=0){
        for (k1 in 1:NROW(out)){
          out$unique.anno.type[k1] <- paste0(sapply(strsplit(out$type[k1],","),unique), collapse=",")
        }
      }
      # count the DMR overlaps; the output can be used to make a barplot or pie-chart
      # unique annotation overlaps

      out.1 <- out[which(sapply(strsplit(out$type,','), uniqueN)==1),]
      out.2 <- out[which(sapply(strsplit(out$type,','), uniqueN)>=2),]

      #counting unique annotations
      for (k2 in seq_along(annotation)){
        anno.list[[k2]] <- NROW(out.1[grep(annotation[k2], out.1$type),])
        names(anno.list)[[k2]] <- annotation[k2]
      }
      #also counting multiple overlaps
      anno.list[length(annotation)+1] <- list(NROW(out.2))
      names(anno.list)[[length(annotation)+1]] <- "multiple.overlaps"
      df.1 <- do.call(cbind, anno.list)

      df <- cbind(sample=tmp.name, total.DMRs=NROW(gr), df.1)
      final.df <- rbind(final.df, data.frame(df))

      out <- out[,-c(4,6)]
      fwrite(x=out, file=paste0(out.dir, "/", tmp.name, "_annotation.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    }
    fwrite(x=final.df, file=paste0(out.dir, "/DMR-counts.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  }
}

