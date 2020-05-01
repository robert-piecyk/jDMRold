makeDMRmatrix<-function(sampleKey, context, chr)
{
  #context="CG"
  #chr=1
  #sampleKey = sampleKey
  all.files<-list.files("DATA/PRODUCED/")
  status.all<-NULL
  level.all<-NULL
  loc.info.all<-NULL
  for (i0 in 1:length(chr))
  {
    chr.temp<-chr[i0]
        status.collect<-NULL
        level.collect<-NULL
            for (i1 in 1:nrow(sampleKey))
            {
              temp<-all.files[grep(as.character(sampleKey[i1, "Code"]), all.files)]
              temp<-temp[grep(paste0("chr", chr.temp, "_", context), temp)]
              temp<-data.frame(fread(paste0("DATA/PRODUCED/", temp)))
              status.collect<-cbind(status.collect, ifelse(as.character(temp[,"status"]) == "U", 0, 1))
              level.collect<-cbind(level.collect, temp[,"rc.meth.lvl"])
            }
            level.collect<-as.matrix(data.frame(level.collect))
            colnames(status.collect)<-sampleKey[,2]
            colnames(level.collect)<-sampleKey[,2]
            loc.info<-temp[,c("seqnames", "start", "end", "context")]
            index<-which(rowSums(status.collect) != 0 & rowSums(status.collect) != ncol(status.collect))
            status.collect<-status.collect[index,]
            level.collect<-level.collect[index,]
            loc.info<-loc.info[index,]
            #index<-lengths(apply(status.collect, 1, table))
            #status.collect<-status.collect[which(index > 1),]
            #level.collect<-level.collect[which(index > 1),]
           # loc.info<-loc.info[which(index > 1),]
      status.all<-rbind(status.all, status.collect)
      level.all<-rbind(level.all, level.collect)
      loc.info.all<-rbind(loc.info.all, loc.info)
  }
  out<-list(status.all, level.all, loc.info.all)
  names(out)<-c("status", "level", "loc.info")
  out
}
filterDMRmatrix<-function(datain, treatment, control)
{
  #treatment<-1:2
  #control<-7:8
  #datain<-df
  level.in<-datain[[2]][,c(treatment, control)]
  status.in<-datain[[1]][,c(treatment, control)]
  loc.in<-datain[[3]]
  index<-which(rowSums(status.in) != 0 & rowSums(status.in) != ncol(status.in))
  status.in<-status.in[index,]
  level.in<-level.in[index,]
  loc.in<-loc.in[index,]
  treat.status<-status.in[,1:length(treatment)]
  control.status<-status.in[, c(length(treatment)+1):ncol(status.in)]
  index1<-which(rowSums(treat.status) == 0 | rowSums(treat.status) == ncol(treat.status))
  index2<-which(rowSums(control.status) == 0 | rowSums(control.status) == ncol(control.status))
  index<-sort(intersect(index1, index2))
  status<-status.in[index,]
  loc<-loc.in[index,]
  level<-level.in[index,]
  out<-list(status, level, loc)
  names(out)<-c("status", "level", "loc.info")
  out
}