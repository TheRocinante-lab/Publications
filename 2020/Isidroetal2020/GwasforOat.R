

setwd("~/Dropbox/AdventuresInMultiOmicsTwo")
load("~/Dropbox/AdventuresInMultiOmicsTwo/Datareadybefore_Gwas.RData")
source('~/Dropbox/AdventuresInMultiOmicsTwo/distancebasedconcordance.R', echo=TRUE)
source('~/Dropbox/AdventuresInMultiOmicsTwo/kernelbasedconcordance.R', echo=TRUE)
source('~/Dropbox/AdventuresInMultiOmicsTwo/HaploFind.R', echo=TRUE)
source('~/Dropbox/AdventuresInMultiOmicsTwo/HaploTest.R', echo=TRUE)
source('~/Dropbox/AdventuresInMultiOmicsTwo/TestTree.R')


#Pheno<-read.csv("BLUPS_gwasoats_data.csv", row.names=1)
Pheno <- pheno
genoimputed[1:5,1:5]
Pheno[1:5,1:11]


#BLUPs
library(lme4)
str(Pheno)

fit.blups <- lmer(HT2_Trans~Year+Block+ (1|GID2), data=Pheno)
BLUPS_HT <- ranef(fit.blups)
Pheno <- data.frame(GID=rownames(BLUPS_HT$GID2),BLUPS=BLUPS_HT$GID2$`(Intercept)`)

traitname <- "Height"
fitBLUP <- function(traitname){
  formula <- paste(traitname,"~ (1|GID2)",sep="")
  formula <- as.formula(formula)
  fit.blups <- lmer(formula, data=pheno)
  BLUPS_trait <- ranef(fit.blups)
  Pheno <- data.frame(GID=rownames(BLUPS_trait$GID2),BLUPS=BLUPS_trait$GID2$`(Intercept)`)
  Pheno$traitname <- traitname
  return(Pheno)
}
colnames(pheno) <- make.names(colnames(pheno))
str(pheno)
BLUPS <- lapply(colnames(pheno)[c(12:24,26:39)], fitBLUP)

library(reshape2)
library(tidyverse)
BLUPS <- Reduce("rbind",BLUPS)

BLUPS <- BLUPS %>% 
pivot_wider(names_from = traitname,values_from=BLUPS)
Pheno <- as.data.frame(BLUPS)

Pheno<-Pheno[match(rownames(genoimputed),Pheno$GID),]
head(map)
colnames(map)<-c("rs","chr","pos")


Map<-map
table(Map$chr)
Map$chr<-as.numeric(gsub("Mrg","",Map$chr))
unique(Map$chr)


genoimputed<-genoimputed[,match(Map$rs,colnames(genoimputed))]

chrorder<-sort(unique(unique(Map$chr)))
DendList<-lapply(sort(unique(unique(Map$chr))),function(x){
  datadist<-data.frame(pos=Map[Map$chr==x,]$pos)
  rownames(datadist)<-Map[Map$chr==x,]$rs
  as.dendrogram(hclust(dist(datadist)))})
library(data.tree)


depth(DendList[[1]])

plot(DendList[[1]])

length(chrorder)

ALLGenoDend<-merge(DendList[[1]],DendList[[2]],DendList[[3]],DendList[[4]],DendList[[5]]
                   ,DendList[[6]],DendList[[7]],DendList[[8]],DendList[[9]],DendList[[10]]
                   ,DendList[[11]],DendList[[12]],DendList[[13]],DendList[[14]],DendList[[15]]
                   ,DendList[[16]],DendList[[17]],DendList[[18]],DendList[[19]],DendList[[20]]
                   ,DendList[[21]])

Pheno <- Pheno[,apply(Pheno,2,function(x){sum(x==0)==0})]
colnames(Pheno)
outSum_t2_HT2<-TestTree(ALLGenoDend, genoimputed,as.data.frame(Pheno[,c(17)]), nlevels = 100, alpha=.05/5000)
out_T2<-TestTree(ALLGenoDend, genoimputed,as.data.frame(Pheno[,c(18)]), nlevels = 100, alpha=.05/5000)
out_H2<-TestTree(ALLGenoDend, genoimputed,as.data.frame(Pheno[,c(19)]), nlevels = 100, alpha=.05/5000)
out_height<-TestTree(ALLGenoDend, genoimputed,as.data.frame(Pheno[,c(2)]), nlevels = 100, alpha=.05/5000)

dfTT<-data.tree::ToDataFrameTree(outSum_t2_HT2, "level", "pval", "isLeaf")
dim(dfTT)
pval<-paste("c(",dfTT$pval,")",sep="")
pval[[1]]
pval<-lapply(1:nrow(dfTT), function(i)eval(parse(text=pval[[i]])))
pval<-Reduce("rbind",pval)
dim(pval)

pval[1:5,1:5]
dfTT<-data.tree::ToDataFrameTree(out, "level", "pval", "isLeaf")


library(tidyverse)
library(ggrepel)
source("manhattanfunction.R")
colnames(Pheno)

plot.fun <- function(x,thershold, ylim,traitname){
  
  dfTT<-data.tree::ToDataFrameTree(x, "level", "pval", "isLeaf")
  pval<-paste("c(",dfTT$pval,")",sep="")
  pval[[1]]
  pval<-lapply(1:nrow(dfTT), function(i)eval(parse(text=pval[[i]])))
  pval<-Reduce("rbind",pval)
  dfTT<-data.tree::ToDataFrameTree(x, "level", "pval", "isLeaf")
  dfTT<-cbind(dfTT,pval)
  dfTT<-na.omit(dfTT)
  AllSNPs <- dfTT[dfTT$isLeaf,]
  AllSNPs$SNPs <- gsub("[[:space:]]","",AllSNPs$levelName)
  AllSNPs$SNPs <- gsub("¦","",AllSNPs$SNPs)
  AllSNPs$SNPs <- gsub("-","",AllSNPs$SNPs)
  AllSNPs$SNPs <- gsub("°","",AllSNPs$SNPs)

  AllSNPs <- AllSNPs[match(Map$rs,AllSNPs$SNPs),]
  Map$pval <- -log10(AllSNPs$pval)
  plots <- manhattan_fun(Map,thershold,ylim,traitname)
  return(plots)
  
}


plot.fun(outSum_t2_HT2,4.29,c(0,40),"Sum_transf")
plot.fun(out_T2,4.29,c(0,40),"T2_Transf")
plot.fun(out_H2,4.29,c(0,40),"HT2_Trans")
plot.fun(out_height,4.29,c(0,40),"Height")


