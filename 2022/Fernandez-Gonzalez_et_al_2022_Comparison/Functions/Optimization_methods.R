################################################################################
################################***CDMEAN***####################################
################################################################################

#For untargeted optimization target population 
#Make sure that the remaining set is always bigger than 50. 
#Otherwise, it is better to calculate CDmean over all candidates
CDmean_without_projection<-function(soln, Data){
  K<-Data[["K"]]
  lambda<-Data[["lambda"]]
  targ <- Data[['Target']]
  if (is.null(targ)) { #for untargeted optimization leave "Target" empty in the data
    cand <- Data[["cand"]]
    targ <- setdiff(cand, soln)
  }
  targ_soln <- c(targ, soln)
  Vinv<-solve(K[soln,soln]+lambda*diag(length(soln))) #inverse of the covariance matrix in the training.
  outmat<-(K[targ_soln,soln]%*%(Vinv-(Vinv%*%Vinv)/sum(Vinv))%*%K[soln,targ_soln])/K[targ_soln,targ_soln] #CD matrix
  
  CD <- mean(diag(as.matrix(outmat[1:length(targ),1:length(targ)]))) #target = test
  
  return(CD)
}



################################################################################
##############################***PCA_CDMEAN***##################################
################################################################################


#used for dimensionality reduction


library(STPGA)

PCA_CDmean<-function(soln, Data){
  PCs_used <- Data[["PCs_used"]]
  data_PCs <- Data[["PCs"]][,1:PCs_used]
  Target <- Data[["Target"]]
  
  PCA_CD <- CDMEAN2(Train = rownames(data_PCs)[soln], Test = rownames(data_PCs)[Target], 
          P = data_PCs, lambda = 1/PCs_used, C=NULL)

  #the minus sign is used because PCA_CD has to be minimized and TrainSel maximizes it
  return(-PCA_CD)
}



################################################################################
##################################***Rscore***##################################
################################################################################

library(Rcpp)

#"Target" is the candidate set in untargetted optimization 
#"Target" is test set in targetted optimization
Rcpp::sourceCpp('./R_Score.cpp')
#dataRopt<-list(X=PC)
optR<-function(soln, Data){
  Xtrain<-Data[["X"]][soln,]
  Xtest<-Data[["X"]][Data[["Target"]],]
  r_score(Xtrain, Xtest)
}


################################################################################
##################***#Average genomic relationship***###########################
################################################################################



#K is the relationship matrix
#"Target" is the candidate set in untargeted optimization 
#"Target" is test set in targeted optimization
Avg_GRM_function<-function(soln, Data){
  K<-Data[["K"]]
  Ksoln<- K[soln,Data[["Target"]]]
  KsolnVec1<- mean(c(unlist(Ksoln)))
  return(c(KsolnVec1))
}

#Same as above but apart from maximizing the relationship TRS-target the 
#average relationship within the TRS is also minimized
Avg_GRM_MinMax_function<-function(soln, Data){
  K<-Data[["K"]]
  target <- Data[["Target"]]
  target <- target[which(!(target %in% soln))] #remove the individuals from the TRS in the target set
  Ksoln <- K[soln,soln]
  Ktarget <- K[soln,target]
  #maximize the relatinship between the TRS (soln) and the target population
  #and minimize the relationship between different elements of the TRS (soln)
  KsolnVec1<- mean(c(unlist(Ktarget))) - mean(c(unlist(Ksoln))) 
  #We could weight differently mean(c(unlist(Ktarget))) and mean(c(unlist(Ksoln))) 
  return(c(KsolnVec1))
}


#Minimize average relationship between individuals in the TRS
#only possible for untargeted optimization
Avg_GRM_self_function<-function(soln, Data){
  K<-Data[["K"]]
  Ksoln<- K[soln,soln]
  #use a negative sign before mean 
  #That way, the average relationship between individuals in the TRS (soln) will be minimized
  KsolnVec1<- -mean(c(unlist(Ksoln)))
  return(c(KsolnVec1))
}



