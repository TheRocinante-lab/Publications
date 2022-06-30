
library(foreach)
library(doMC)
library(rrBLUP)
#library(TrainSel)
#library(emmeans)
library(BGLR)

#Optimization_method_format = "Standard" --> RAND, PAM. This is the default format
#Optimization_method_format = "TrainSel" --> CD, PEV, U and its targetted versions

#pheno must have 2 columns, GID and Trait

Models <- function(geno, 
                   K, 
                   pheno = NULL, 
                   Optimization_method_name, 
                   Optimization_method_format, 
                   selected_train, 
                   selected_test, 
                   niter, 
                   training_set_sizes, 
                   nenvironments, 
                   nreps, 
                   Bayes_burnIn,
                   Bayes_nIter,
                   data_name,
                   trait_name,
                   n.cpus = 1,
                   simulated_trait = NULL,
                   Models_used) {
  

  registerDoMC(cores = n.cpus) 
  
  
  GBLUP <- c()
  BayesB <- c()
  RKHS <- c()
  models <- list(GBLUP = c(), BayesB = c(), RKHS = c())

  #check that there is only a real trait or a simulated trait but not both
  if (!is.null(pheno) && !is.null(simulated_trait)) {
    stop("The input can be a real or a simulated trait, but not both\n")
  }

  
  #need package foreach
  foreach(i = 1:niter)%dopar%{
    
    #load simulated trait if needed
    if (!is.null(simulated_trait)) {
      pheno <- simulated_trait[[i]][["pheno"]]
    }
    
 
        #save the individuals selected as test set in the rep number i.
    test <-rownames(K)[selected_test[[i]]] 


    for(ntr in 1:training_set_sizes) {
      
      cat("\n####################\n")
      cat("Optimization method:", Optimization_method_name, "\n")
      cat("Iteration number:", i, "\n")
      cat("Training set size number:", ntr, "\n")
      cat("####################\n\n")
      
      

      if (Optimization_method_format == "TrainSel") {
        train <- rownames(K)[selected_train[[i]][[ntr]][1]$BestSol_int]
      } else if (Optimization_method_format == "Entire_CS") {
        train <- rownames(K)[selected_train[[i]][[1]]]
      } else {
        train <- rownames(K)[selected_train[[i]][[ntr]]]
      }

      pheno.train <- pheno[pheno$GID%in%train,]
      pheno.test <- pheno[pheno$GID%in%test,]
      geno.train <- geno[rownames(geno)%in%train,]
      geno.test <- geno[rownames(geno)%in%test,]

      train_plus_test <- pheno$GID%in%c(train,test)

     # Xtrain <- model.matrix(~pheno.train$ENV); dim(Xtrain)
    #  Xtest <- model.matrix(~pheno.test$ENV); dim(Xtest)
      
      #I need the X matrix for RHKS, but as I only have BLUPs and there is no
      #information about the environments, the X matrix will be a column full of 1s
     # Xtrain<-matrix(1, nrow=length(train),ncol=1); dim(Xtrain)
     # Xtest<-matrix(1, nrow=length(test),ncol=1); dim(Xtest)
      

      Ztrain <- model.matrix(~pheno.train$GID-1); dim(Ztrain)
      colnames(Ztrain) <- colnames(K)
      Ztrain_train_plus_test <- Ztrain[,train_plus_test]; dim(Ztrain_train_plus_test)
      
      Ztest <- model.matrix(~pheno.test$GID-1); dim(Ztest)
      colnames(Ztest) <- colnames(K)
      Ztest_train_plus_test <- Ztest[,train_plus_test]; dim(Ztest_train_plus_test)
      
      ##########################################################################
      ################################GBLUP#####################################
      ##########################################################################
      
      if ("GBLUP" %in% Models_used) {
      
      cat("\n GBLUP \n")
      
      #train model
      
      #Include only the individuals in the training and test
      #sets to reduce dimensionality and accelerate the calculation
      
      twostep_gblup <- mixed.solve(pheno.train$Trait,
                                   X=NULL,
                                   Z=Ztrain_train_plus_test,
                                   K=K[train_plus_test,train_plus_test])
      
      GEBVs_GBLUP <- twostep_gblup$u
      
      Ve_GBLUP <- twostep_gblup$Ve
      Vu_GBLUP <- twostep_gblup$Vu
      
      
      #Narrow sense
     # h2_GBLUP <- Vu_GBLUP/(Vu_GBLUP+(Ve_GBLUP/(nenvironments*nreps)))
      
      #accuracy
      #cor(Ztest%*%GEBVs_GBLUP,pheno.test$Trait)
      accuracy_GBLUP <- cor(Ztest_train_plus_test%*%GEBVs_GBLUP,pheno.test$Trait)
      
      
      GBLUP[[ntr]] <- list(model = twostep_gblup, 
                          # h2 = h2_GBLUP, 
                           accuracy = accuracy_GBLUP)
      
      models$GBLUP <- GBLUP
      
      
      #free memory
      rm(twostep_gblup)
      #gc = garbage collector --> frees memory from removed objects
      gc()
      
      }
      
      ##########################################################################
      ################################BayesB####################################
      ##########################################################################
      
      if ("BayesB" %in% Models_used) {
      
      cat("\n BayesB \n")
      
      
      # length(pheno.train$Trait)
      # dim(Xtrain)
      
      twostep_BayesB<-BGLR(pheno.train$Trait,
                           ETA=list(list(X=geno.train,model="BayesB")),nIter=Bayes_nIter, 
                           burnIn=Bayes_burnIn,verbose=FALSE,saveAt = paste("Iter_",i,"_",sep=""))
      #In saveAt a different name for each iteration was chosen
      #It is important because otherwise all iterations use the same files for BGLR
      


      
      #ETA is the functions to fit the model. You add design matrix and the
      #type of effect. Fixed or Random. If it´s random which Bayes it is.
      #nIter from the posterior (start from the prior and then each iteration goes close to the posterior)
      #burnIN delete the first 5000 iterations because we assume we are not in the posterior in the first
      #iterations. 
      #verbose just if you want to see the iterations (print) on the screen. 
      
      
      #I do not save marker.effects.BayesB to save memory 
      #marker.effects.BayesB<- twostep_BayesB$ETA[[1]]$b
      
      
      
      # plot(abs(marker.effects.BayesB), ylab='Estimated Squared-Marker Effect',
      #      type='o',cex=.5,main='Marker Effects')
      
      ########After you get marker effects we can get GEBVs
      
      # For the train set lines
      GEBVstrain.BayesB<-(Ztrain%*%geno)%*%twostep_BayesB$ETA[[1]]$b
      ## For the test set lines
      GEBVstest.BayesB<-(Ztest%*%geno)%*%twostep_BayesB$ETA[[1]]$b
      ### For all the lines
      #GEBVs.BayesB<-as.matrix(genodata)%*%marker.effects.BayesB
      
      
      
      ########
      accuracy_BayesB <- cor(pheno.test$Trait, GEBVstest.BayesB)
      
      
      BayesB[[ntr]] <- list(model = twostep_BayesB, 
                            #marker_effects = marker.effects.BayesB, 
                            accuracy = accuracy_BayesB)
      
      models$BayesB <- BayesB
      

      
      #free memory
      rm(twostep_BayesB)
      #gc = garbage collector --> frees memory from removed objects
      gc()
      
      }
      
      ##########################################################################
      ################################RKHS######################################
      ##########################################################################
      
      if ("RKHS" %in% Models_used) {
      
      cat("\n RKHS \n")
      
      #First, test several values of h and select the best one
      
      #Include in Eucl.distance only the individuals in the training and test
      #sets to reduce dimensionality and accelerate the calculation
      Eucl.distance<-as.matrix(dist(geno[train_plus_test,]))
      h<-1/ncol(geno)
      Amat.RKHS<-exp(-h*(Eucl.distance)^2/(2)) #gaussian kernel
      
      #h is a parameter controlling the rate of decay in the kernel matrix (Amat.RKHS)
      #test several values and select the best
      nmarkers<-ncol(geno)
      hvec<-rep(1/nmarkers,5)*c(1/3,1/2,1,2,3)
      hvec
      RHKS.hfun<-function(h){
        Amat.RKHS<-exp(-h*(Eucl.distance)^2/(2))
        model.RKHS<-mixed.solve(y= pheno.train$Trait, X=NULL, Z=Ztrain_train_plus_test, K=Amat.RKHS)
        return(model.RKHS$LL)
      }
      

      
      
      RKHS.models<-sapply(hvec,RHKS.hfun,simplify=TRUE)
      RKHS.models
      h.best<-hvec[which.max(RKHS.models)]
      
      #create Amat with the best h and fit the model
      Amat.RKHS.best.h<-exp(-h.best*(Eucl.distance)^2/(2))
      
      
      model.RKHS.best.h.twostep<-mixed.solve(y= pheno.train$Trait, 
                                             X=NULL, 
                                             Z=Ztrain_train_plus_test,
                                             K=Amat.RKHS.best.h)
      
      
      
      
      
      
      GEBVstrain.RKHS.best.h.twostep<-Ztrain_train_plus_test%*%model.RKHS.best.h.twostep$u
      
      
      GEBVstest.RKHS.best.h.twostep<-Ztest_train_plus_test%*%model.RKHS.best.h.twostep$u
      
      accuracy_RKHS <- cor(pheno.test$Trait, GEBVstest.RKHS.best.h.twostep)
      
      
      
      RKHS[[ntr]] <- list(model = model.RKHS.best.h.twostep, 
                          h_best = h.best, 
                          accuracy = accuracy_RKHS)
      
      models$RKHS <- RKHS
      
      
      #free memory
      rm(RKHS.models)
      rm(model.RKHS.best.h.twostep)
      #gc = garbage collector --> frees memory from removed objects
      gc()
      
      }
   
      ##########################################################################
      ##########################################################################
      ##########################################################################
      
      
      
      #do not iterate for every TRS size when making models for the entire 
      #candidate set as in that case the TRS size is irrelevant
      if (Optimization_method_format == "Entire_CS") {
        break
      }
      
      
    }
    
    if (is.null(simulated_trait))  {
      save(models, file = paste(data_name, "_trait_", trait_name, "_models_", Optimization_method_name, "_rep_", i, ".RData", sep= ""))
    } else {
      #if the trait is simulated, save the maker effects and the phenotype
      heritability_narrow <- simulated_trait[[i]][["h2"]]
      save(models, heritability_narrow, file = paste(data_name, "_trait_", trait_name, "_models_", Optimization_method_name, "_rep_", i, ".RData", sep= ""))
    }
  }

}

