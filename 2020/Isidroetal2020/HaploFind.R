# #######################################
#  
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#  
# BiocManager::install("Biobase")
# library(devtools)
# install_github("sunnyeesl/BigLD")
#  library(BigLD)
#  data(geno)
#  data(SNPinfo)


############



HaploFind<-function(Geno, map.geno, ...){
  intersectnames<-intersect(colnames(Geno), map.geno$rs)
  map.geno<-map.geno[map.geno$rs%in%intersectnames,]
  Geno<-Geno[, colnames(Geno)%in%intersectnames]
  Geno<-Geno[,match(map.geno$rs,colnames(Geno))]
  ###geno is a dataframe with columns as snps
  ###  snp info is dataframe with colums rsID, bp
  MAPSPLITCHR<- split.data.frame(map.geno, as.factor(as.character(map.geno$chr)))
  GenoSPLITCHR<- split.data.frame(t(Geno),  as.factor(as.character(map.geno$chr)))
  
  GenoSPLITCHR<-lapply(GenoSPLITCHR, function(x){as.data.frame(t(x))})
  outlapply<-lapply(1:length(MAPSPLITCHR), function(i){
    genochr<-GenoSPLITCHR[[i]]
    SNPinf<-MAPSPLITCHR[[i]]
    SNPinf<-SNPinf[,colnames(SNPinf)%in%c("rs","pos")]
    SNPinf<-SNPinf[,c("rs","pos")]
    colnames(SNPinf)<-c("rsID", "bp")
    haploout<-BigLD::Big_LD(genochr, SNPinf, ...)
    poshaploout<-c(unlist(apply(haploout,1,function(x){c(x[1]:x[2])})))
    remainingrs<-setdiff(1:nrow(SNPinf),poshaploout)
    remainingmap<-SNPinf[remainingrs,]
    remainingmap$row<-remainingrs
    remaining<-Reduce("rbind", lapply(remainingmap$row,function(i){
      return(as.data.frame(cbind(start=i, end=i,  start.rsID=as.character(remainingmap[remainingmap$row==i,]$rsID),   end.rsID=as.character(remainingmap[remainingmap$row==i,]$rsID), start.bp=remainingmap[remainingmap$row==i,]$bp,   end.bp=remainingmap[remainingmap$row==i,]$bp)))
    }))
    haploout<-rbind(haploout,remaining)
    haploout<-haploout[order(as.numeric(haploout$start)),]
    haploout$haplo<-1:nrow(haploout)
    rownames(haploout)<-NULL
    return(haploout)
  }
  )
  
  return(outlapply)
}



Haplo2Names<-function(hapout,map.geno, Geno){
  
  intersectnames<-intersect(colnames(Geno), map.geno$rs)
  map.geno<-map.geno[map.geno$rs%in%intersectnames,]
  ###geno is a dataframe with columns as snps
  ###  snp info is dataframe with colums rsID, bp
  MAPSPLITCHR<- split.data.frame(map.geno, as.factor(as.character(map.geno$chr)))
  hapout<-Reduce("rbind",lapply(1:length(MAPSPLITCHR), function(i){
    print(i)
    mapi<-MAPSPLITCHR[[i]]
    mapi<-mapi[mapi$rs%in%intersectnames,]
    hapvec<-c()
    for (j in 1:nrow(hapout[[i]])){
      hapvec<-c(hapvec,rep(hapout[[i]][j,]$haplo, as.numeric(hapout[[i]][j,2])-as.numeric(hapout[[i]][j,1])+1))
    }
    mapi$haplo<-hapvec
    return(mapi)
  }))
  
  return(hapout)
}



meanposvec<-function(map, annotation){
  aggout<-aggregate(map$pos, by=list(chr=map$chr, ann=annotation), mean)
  x<-aggout$x[match(paste(map$chr,annotation),paste(aggout$chr, aggout$ann))]
  return(x)
}



HaploFindPhyPos<-function(Geno, map.geno, numparts=10){
  intersectnames<-intersect(colnames(Geno), map.geno$rs)
  map.geno<-map.geno[map.geno$rs%in%intersectnames,]
  Geno<-Geno[, colnames(Geno)%in%intersectnames]
  Geno<-Geno[,match(map.geno$rs,colnames(Geno))]
  ###geno is a dataframe with columns as snps
  ###  snp info is dataframe with colums rsID, bp
  MAPSPLITCHR<- split.data.frame(map.geno, as.factor(as.character(map.geno$chr)))
  GenoSPLITCHR<- split.data.frame(t(Geno),  as.factor(as.character(map.geno$chr)))
  
  GenoSPLITCHR<-lapply(GenoSPLITCHR, function(x){as.data.frame(t(x))})
  outlapply<-Reduce("rbind",lapply(1:length(MAPSPLITCHR), function(i){
    genochr<-GenoSPLITCHR[[i]]
    SNPinf<-MAPSPLITCHR[[i]]
    SNPinf<-SNPinf[,colnames(SNPinf)%in%c("rs","chr","pos")]
    SNPinf<-SNPinf[,c("rs","chr","pos")]
    colnames(SNPinf)<-c("rsID", "chr","bp")
    SNPinfsplit<-split(SNPinf, cut(SNPinf$bp, numparts))
    SNPinfsplit<-Reduce("rbind",lapply(1:length(SNPinfsplit), function(j){
      xx<-SNPinfsplit[[j]]
      if (nrow(xx)>0){
      xx$haplo<-j
      } else {
        xx$haplo<-c()
      }
      return(xx)
      }))
    return(SNPinfsplit)
  }
  ))
  
  return(outlapply)
}





HaploFindPhyPosHyClust<-function(Geno, map.geno, numparts=10){
  intersectnames<-intersect(colnames(Geno), map.geno$rs)
  map.geno<-map.geno[map.geno$rs%in%intersectnames,]
  Geno<-Geno[, colnames(Geno)%in%intersectnames]
  Geno<-Geno[,match(map.geno$rs,colnames(Geno))]
  ###geno is a dataframe with columns as snps
  ###  snp info is dataframe with colums rsID, bp
  if (length(numparts)==1){numparts<-rep(numparts,length(unique(map.geno$chr)))}
  MAPSPLITCHR<- split.data.frame(map.geno, as.factor(as.character(map.geno$chr)))
  GenoSPLITCHR<- split.data.frame(t(Geno),  as.factor(as.character(map.geno$chr)))
  GenoSPLITCHR<-lapply(GenoSPLITCHR, function(x){as.data.frame(t(x))})
  outlapply<-Reduce("rbind",lapply(1:length(MAPSPLITCHR), function(i){
    genochr<-GenoSPLITCHR[[i]]
    SNPinf<-MAPSPLITCHR[[i]]
    SNPinf<-SNPinf[,colnames(SNPinf)%in%c("rs","chr","pos")]
    SNPinf<-SNPinf[,c("rs","chr","pos")]
    colnames(SNPinf)<-c("rsID", "chr","bp")
    Dchr<-dist(SNPinf$bp)
    hcl<-hclust(Dchr)
    library(ape)
    clustering<-cutree(hcl, k=numparts[i])
    #plot(hcl)
    
    plot(as.phylo(hcl), type = "unrooted", cex = 0.6,
         no.margin = TRUE, tip.color=clustering, main=paste("chr",i))
    #points(.5,.5,col="red")
    SNPinf$haplo<-clustering
    return(SNPinf)
  }
  ))
  return(outlapply)
}







  