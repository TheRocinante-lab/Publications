setwd("~/Dropbox/GenomicMating/GenomicMatingArticleV2/GM_v2.0")
rm(list=ls())
load("wheatdataready.RData") ## data from training population CHKS data, the pheno are blups already.


##### The idea is to first do Training population optimization and then GM.
library(STPGA)

niterations=3
nelite=5
npop=10
mutprob=.8
mutintensity=2
trainsizevec=c(25)
errorstatvec<-c("CDMEAN","DOPT","PEVMEANMM")
errorstat=errorstatvec[1]


ListTrain1<-GenAlgForSubsetSelectionNoTest(P=geno,ntoselect=trainsizevec,npop=npop, nelite=nelite, 
                                           mutprob=mutprob, mutintensity=mutintensity,
                                           niterations=niterations, lambda=1e-5, 
                                           plotiters = T, errorstat=errorstat,C=NULL, 
                                           mc.cores=3, tolconv=1e-5)

ListTrain2<-GenAlgForSubsetSelection(P=PCgeno,
                                     Candidates=rownames(PCgeno)[!(rownames(PCgeno)%in%LeaveoutforTest)], Test=LeaveoutforTest, ntoselect=ntrain,npop=npop, nelite=nelite, mutprob=mutprob, mutintensity=mutintensity,niterations=niterations, lambda=1e-5, plotiters = T, errorstat=errorstat,C=NULL, mc.cores=3, tolconv=1e-5)


allselectected[[itrainselect]]<-ListTrain1
itrainselect=itrainselect+1
allselectected[[itrainselect]]<-ListTrain2
itrainselect=itrainselect+1
minimums<-c(minimums,min(ListTrain1[[nelite+1]]))
minimums2<-c(minimums2,min(ListTrain2[[nelite+1]]))
train<-factor(ListTrain1[[1]], levels=rownames(amatout$A))
train2<-factor(ListTrain2[[1]], levels=rownames(amatout$A))
test=factor(setdiff(rownames(X),train), levels=rownames(amatout$A))
test2=factor(setdiff(rownames(X),train2), levels=rownames(amatout$A))

LeaveoutforTest=factor(LeaveoutforTest, levels=rownames(amatout$A))

Ztrain<-model.matrix(~-1+train)
Ztrain2<-model.matrix(~-1+train2)
Ztest<-model.matrix(~-1+test)
Ztest2<-model.matrix(~-1+test2)
ZtestLeaveout<-model.matrix(~-1+LeaveoutforTest)





gasols<-getGaSolutions(Markers=geno,Markers2=NULL, K=A.mat, markereffects=markereffects,markermap=NULL,nmates=10,
                       minparents=3, 
                       impinbreedstepsize=.02, impvar=.01, 
                       impforinbreed=.01,
                       npopGA=100, nitGA=10, miniters=10,minitbefstop=20,plotiters=T,
                       mc.cores=1,nelite=20, mutprob=0.8, noself=T,
                       method=1, type=0L, generation=0L)


