HaploTest<-function(mapwithhaplo, Pheno, Geno){
  pdout<-NULL 
  x=dist(Pheno)

for (chr in sort(as.numeric(unique(mapwithhaplo$chr)))){
 print(chr)
  mapwithhaplochr<-mapwithhaplo[mapwithhaplo$chr==chr,]
  for (haplo in unique(mapwithhaplochr$haplo)){
    print(haplo)
  inmarkers<-which((mapwithhaplo$chr==chr & mapwithhaplo$haplo==haplo))
  if (length(inmarkers)>0){
  outmarkers<-setdiff(1:ncol(Geno), inmarkers)
  y=Geno[,inmarkers]
  if (length(outmarkers)>0){
  z=Geno[,outmarkers]
  pdout<-rbind(pdout,partialdistcov(x,y,z,2)$pdcortest)
  } else {
    pdout<-rbind(pdout,distcov(x,y)$dcortest)
  }
  } else {
    pdout<-rbind(pdout,1)
  }
  }
}
return(pdout)
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

