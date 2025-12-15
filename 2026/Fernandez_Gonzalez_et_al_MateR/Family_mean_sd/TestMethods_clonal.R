set.seed(1234)
n<- 0 
n2 <- 10000
nrep <- 500
DH=FALSE
phi=2
nQTL = 250


library(MateR)
library(SimpleMating)
library(predCrossVar)
data("ExampleDataDiploid")

# introduce your MateR license keys here.
# To obtain license keys, please contact j.isidro@upm.es.
# Free license keys will be provided to public bodies,
# such as, Universities, and non-profit organizations.
Username=username
Password=password


#This value will be used for all chromosomes.
ChromosomeCentiMorgan <- 150
ChromosomesCentiMorgan <- rep(ChromosomeCentiMorgan, length(Chromosomes))

generic_GenMap_SM <- c()
for (i in 1:length(Chromosomes)) {
  size <- ChromosomesCentiMorgan[[i]]
  nMarkers <- length(Chromosomes[[i]])
  interval <- size/nMarkers
  
  tmp <- data.frame(chr = rep(i, nMarkers),
                    pos = seq(interval, size, by=interval),
                    mkr = names(markereffects$YLD[Chromosomes[[i]]]))
  generic_GenMap_SM <- rbind(generic_GenMap_SM, tmp)
}


Chromosomes_names <- Chromosomes
Chromosomes<-lapply(Chromosomes, function(x){return(which(colnames(Markers)%in%x))})



#Follow Haldane Model (number of recombination events follow Poisson distribution)
create_gamete <- function(haplotype, Chromosomes, phi=2, ChromosomesCentiMorgan) {
  
  if (phi == 2) {
    gamete <- c()
    for (chr in 1:length(Chromosomes)) {
      chromosome <- Chromosomes[[chr]]
      cM <- ChromosomesCentiMorgan[chr]
      order <- sample(1:phi,phi)
      nrecombinations <- rpois(1, lambda = cM/100)
      recombination_positions <- sort(sample(chromosome[-1],nrecombinations), decreasing = F)
      if (nrecombinations == 0) {
        gamete <- c(gamete, haplotype[order[1],chromosome[1]:chromosome[length(chromosome)]])
      } else {
        for (p in 0:nrecombinations) {
          if (p == 0) {
            gamete <- c(gamete, haplotype[order[p%%2+1],chromosome[1]:(recombination_positions[p+1]-1)])
          } else if (p == nrecombinations) {
            gamete <- c(gamete, haplotype[order[p%%2+1],(recombination_positions[p]):chromosome[length(chromosome)]])
          } else {
            gamete <- c(gamete, haplotype[order[p%%2+1],recombination_positions[p]:(recombination_positions[p+1]-1)])
          }
        }
      }
    }
  }
  
  if (phi == 4) {
    #Assume only bivalents in chiasmata, i.e., the 4 chromosomes are organized
    #into 2 pairs that only recombinate within themselves. 
    #This is frequent in autotetraploids: https://link.springer.com/article/10.1007/s00412-015-0571-4
    order <- sample(1:phi,phi)
    bivalent1 <- haplotype[order[1:2],] #first pair of chromosomes
    bivalent2 <- haplotype[order[3:4],] #second pair of chromosomes
    
    #treat the bivalents as a diploid
    gamete1 <- create_gamete(bivalent1, Chromosomes, phi=2, ChromosomesCentiMorgan) 
    gamete2 <- create_gamete(bivalent2, Chromosomes, phi=2, ChromosomesCentiMorgan) 
    gamete <- rbind(gamete1, gamete2)
  }
  
  return(gamete)
}





#adapted from Wolfe GitHub
#https://github.com/wolfemd/predCrossVar/blob/master/R/predCrossVar.R
Wolfe_GVs <- function(sireID,damID,doseMat,postMeanAddEffects,postMeanDomEffects=NULL, p_pop, q_pop){
  if (!is.null(postMeanDomEffects)) {
    
    # Eqn 14.6 from Falconer+MacKay
    p1<-doseMat[sireID,]/2
    p2<-doseMat[damID,]/2
    q<-1-p1
    y<-p1-p2
    g<-postMeanAddEffects*(p1-q-y) + postMeanDomEffects*((2*p1*q)+y*(p1-q))
    meanG<-sum(g)
    # Equivalent, e.g. Toro & Varona 2010
    # q1<-1-p1
    # q2<-1-p2
    # gfreqs<-cbind(q1*q2,p1*q2+p2*q1,p1*p2)
    # meanG<-sum((gfreqs[,3]-gfreqs[,1])*postMeanAddEffects[[Trait]]+gfreqs[,2]*postMeanDomEffects[[Trait]])
    return(meanG)
  } else {
    sireGEBV <- doseMat[sireID,]%*%postMeanAddEffects
    damGEBV <- doseMat[damID,]%*%postMeanAddEffects
    (sireGEBV+damGEBV)/2
  }

}




#GenomicMating var
#Adapted from GenomicMating github and outputting the same values as their C++ code:
getstatsM1R <- function(MARKERS, K, markereffects, P) {
  output <- rep(NA, 3)
  EBV = (MARKERS+1)%*%markereffects #our marker effects are designed for a matrix coded as 0,1,2; but the input of this function is in -1, 0, 1. Undo it so that they match
  output[2] <- rep(1, nrow(P))%*%P%*%EBV #sum of mid-parental values for each cross
  

  output[1] <- NA #this contains inbreeding results that are not needed for
  #within-family computations. We will just skip it.
  
  crossvalues <- list()
  Parents <- matrix(nrow = 2, ncol = ncol(MARKERS))
  
  for (i in 1:nrow(P)) {
    tempbool = which(P[i,] > 0)
    
    if (length(tempbool) == 1) {
      Parents[1,] <- MARKERS[tempbool,]
      Parents[2,] <- MARKERS[tempbool,]
    } else {
      Parents <- MARKERS[tempbool,]
    }
    
    
    absmarkereffects <- abs(markereffects)
    
    p1p2 <- Parents
    
    #this is a much easier way to do everything the package does in c++
    scoresv <- rep(0,ncol(p1p2))
    scoresv[which(p1p2[1,]==0)] <- 0.25
    scoresv[which(p1p2[2,]==0)] <- 0.25
    scoresv[which(p1p2[1,]==0 & p1p2[2,]==0)] <- 0.5
    
    
    ############################################################################
    ############################################################################
    ############################################################################
    #Typo in formula, family variances get multiplied by 4!!!!!!!!!
    #what GenomicMating does
    crossvalues[[i]] <- scoresv*(2*absmarkereffects)^2
    
    #what it should be:
    #crossvalues[[i]] <- scoresv*(absmarkereffects)^2
    ############################################################################
    ############################################################################
    ############################################################################
    
  }
  output[3] <- sum(unlist(crossvalues)) #sum of family variance for each cross
  return(output)
  
}









#Heterozygous parents, F1

additive_values_list <-  list()
genotypic_values_list <-  list()
parents <- list()


#Initialize vectors for results of within-family sd.
#In hindsight, I should have just put them in a list instead of initializing
#so many vectors. I did that for the family mean. 

#naming convention:
#the addition of "hap" indicates that full haplotype data was used to compute LD
#the addition of "sigma" indicates LD was approximated from the parental population
#the addition of "ind" indicates LD was disregarded
#the addition of "dom" indicates dominance effects were considered

#MateR
sd_hap <- c() 
sd_sigma <- c() 
sd_ind <- c() 
sd_hap_dom <- c()
sd_sigma_dom <- c()
sd_ind_dom <- c()

#wolfe 2021 (predCrossVar)
sd_wolfe <- c()
sd_wolfe_dom <- c()

#Simple Mating
sd_hap_SM <- c()
sd_sigma_SM <- c()
sd_hap_dom_SM <- c()
sd_sigma_dom_SM <- c()




#MateR
sd_hap_hat <- c()
sd_sigma_hat <- c()
sd_ind_hat <- c()
sd_hap_dom_hat <- c()
sd_sigma_dom_hat <- c()
sd_ind_dom_hat <- c()

#wolfe 2021 (predCrossVar)
sd_wolfe_hat <- c()
sd_wolfe_dom_hat <- c()

#Simple Mating
sd_hap_SM_hat <- c()
sd_sigma_SM_hat <- c()
sd_hap_dom_SM_hat <- c()
sd_sigma_dom_SM_hat <- c()




#MateR
sd_hap_hat_0.5 <- c()
sd_sigma_hat_0.5 <- c()
sd_ind_hat_0.5 <- c()
sd_hap_dom_hat_0.5 <- c()
sd_sigma_dom_hat_0.5 <- c()
sd_ind_dom_hat_0.5 <- c()

#wolfe 2021 (predCrossVar)
sd_wolfe_hat_0.5 <- c()
sd_wolfe_dom_hat_0.5 <- c()

#Simple Mating
sd_hap_SM_hat_0.5 <- c()
sd_sigma_SM_hat_0.5 <- c()
sd_hap_dom_SM_hat_0.5 <- c()
sd_sigma_dom_SM_hat_0.5 <- c()




#MateR
sd_hap_realistic <- c()
sd_sigma_realistic <- c()
sd_ind_realistic <- c()
sd_hap_dom_realistic <- c()
sd_sigma_dom_realistic <- c()
sd_ind_dom_realistic <- c()

#wolfe 2021 (predCrossVar)
sd_wolfe_realistic <- c()
sd_wolfe_dom_realistic <- c()

#Simple Mating
sd_hap_SM_realistic <- c()
sd_sigma_SM_realistic <- c()
sd_hap_dom_SM_realistic <- c()
sd_sigma_dom_SM_realistic <- c()




#GenomicMating
sd_GM <- c()
sd_GM_hat <- c()
sd_GM_hat_0.5 <- c()
sd_GM_realistic <- c()

means <- list()
for (rep in 1:nrep) {
  
  
  #additive
  QTL_numeric <- sample(1:length(markereffects$YLD), nQTL)
  Marker_positions <- colnames(Markers)[setdiff(1:length(markereffects$YLD), QTL_numeric)]
  QTL <- colnames(Markers)[QTL_numeric]
  
  
  #Get estimated, non-true marker effects
  #additive
  BVs <- c(Markers[,QTL]%*%markereffects$YLD[QTL])
  names(BVs) <- rownames(Markers[,QTL])
  markereffects_hat <- c(meff_from_GVs(Markers[,QTL],
                                       parametrization = "Genotypic",
                                     phi,
                                     BVs,
                                     type = "additive",
                                     DirectionalDominance = NULL))
  names(markereffects_hat) <- names(markereffects$YLD[QTL])
  
  sd_a <- sd(BVs)
  #we want a 0.5 heritability --> accuracy will be sqrt(0.5)
  #sd_a/(sd_a+sd_rest) = 0.5; sd_rest = sd_a
  EBVs <- BVs + rnorm(length(BVs), mean = 0, sd = sd_a)
  #scale back to same sd
  EBVs <- EBVs/sd(EBVs)*sd_a
  sd(EBVs)
  cor(EBVs, BVs)
  
  markereffects_hat_0.5 <- c(meff_from_GVs(Markers[,QTL], 
                                           parametrization = "Genotypic",
                                         phi,
                                         EBVs,
                                         type = "additive",
                                         DirectionalDominance = NULL))
  names(markereffects_hat_0.5) <- names(markereffects$YLD[QTL])
  
  
  cor(markereffects$YLD[QTL], markereffects_hat)
  cor(markereffects$YLD[QTL], markereffects_hat_0.5)
  
  cor(Markers[,QTL]%*%markereffects$YLD[QTL], Markers[,QTL]%*%markereffects_hat)
  cor(Markers[,QTL]%*%markereffects$YLD[QTL], Markers[,QTL]%*%markereffects_hat_0.5)
  
  
  #dominance
  #we need heterozygosity. Create some random F1s
  #No need of recombination here, parental lines are fully homozygous
  Markers_d <- c()
  for (i in 1:100) {
    Markers_d <- rbind(Markers_d, colSums(H_parents[[sample(1:length(H_parents),1)]] +  H_parents[[sample(1:length(H_parents),1)]])/2)
  }
  Markers_d[Markers_d!=1] <- 0
  
  DVs <- c(Markers_d[,QTL]%*%markereffects_d$YLD[QTL])
  names(DVs) <- rownames(Markers_d[,QTL])
  markereffects_d_hat <- c(meff_from_GVs(Markers_d[,QTL], 
                                         parametrization = "Genotypic",
                                       phi,
                                       DVs,
                                       type = "dominance",
                                       DirectionalDominance = NULL))
  names(markereffects_d_hat) <- names(markereffects$YLD[QTL])
  
  
  
  
  sd_d <- sd(DVs)
  #we want a 0.5 heritability --> accuracy will be sqrt(0.5)
  #sd_a/(sd_a+sd_rest) = 0.5; sd_rest = sd_a
  EDVs <- DVs + rnorm(length(DVs), mean = 0, sd = sd_d)
  #scale back to same sd
  EDVs <- EDVs/sd(EDVs)*sd_d
  sd(EDVs)
  cor(EDVs, DVs)
  
  markereffects_d_hat_0.5 <- c(meff_from_GVs(Markers_d[,QTL], 
                                             parametrization = "Genotypic",
                                           phi,
                                           EDVs,
                                           type = "dominance",
                                           DirectionalDominance = NULL))
  names(markereffects_d_hat_0.5) <- names(markereffects$YLD[QTL])
  
  
  cor(markereffects_d$YLD[QTL], markereffects_d_hat)
  cor(markereffects_d$YLD[QTL], markereffects_d_hat_0.5)
  
  cor(Markers_d[,QTL]%*%markereffects_d$YLD[QTL], Markers_d[,QTL]%*%markereffects_d_hat)
  cor(Markers_d[,QTL]%*%markereffects_d$YLD[QTL], Markers_d[,QTL]%*%markereffects_d_hat_0.5)
  
  
  
  
  
  #QTL different than markers
  #realistic scenario
  
  markereffects_realistic <- c(meff_from_GVs(Markers[,Marker_positions], 
                                             parametrization = "Genotypic",
                                           phi,
                                           EBVs,
                                           type = "additive",
                                           DirectionalDominance = NULL))
  names(markereffects_realistic) <- names(markereffects$YLD[Marker_positions])
  
  
  cor(Markers[,QTL]%*%markereffects$YLD[QTL], Markers[,Marker_positions]%*%markereffects_realistic)
  
  
  
  # #dominance

  markereffects_d_realistic <- c(meff_from_GVs(Markers_d[,Marker_positions], 
                                               parametrization = "Genotypic",
                                             phi,
                                             EDVs,
                                             type = "dominance",
                                             DirectionalDominance = NULL))
  names(markereffects_d_realistic) <- names(markereffects$YLD[Marker_positions])
  
  

  cor(Markers_d[,QTL]%*%markereffects_d$YLD[QTL], Markers_d[,Marker_positions]%*%markereffects_d_realistic)
  
  
  
  
  
  
  
  
  
  
  
  
  #Create an F1 without recombination (no need to recombine fully homzygous parents)
  parents1 = "tmp"
  parents2 = "tmp" 
  while (sum(parents1 == parents2)>0) {
    parents1 <- sample(rownames(Markers),2)
    F11 <- rbind(H_parents[[parents1[1]]][1,], H_parents[[parents1[2]]][2,])
    parents2 <- sample(rownames(Markers),2)
    F12 <- rbind(H_parents[[parents2[1]]][1,], H_parents[[parents2[2]]][2,])
  }
  
  parents[[rep]] <- c(parents1, parents2)

  additive_values <- c()
  genotypic_values <- c()
  
  print(rep)
  
  #progress bar
  cat("0%   10   20   30   40   50   60   70   80   90   100%\n", file = stderr())
  cat("[----|----|----|----|----|----|----|----|----|----|\n", file = stderr())
  
  for (i in 1:n2) {
    
    #Progress bar
    counter <- i
    interval <- n2*0.02
    previous <- floor((counter-1)/interval)
    current <- floor(counter/interval)
    if (current > previous) {
      cat(paste(rep("*", (current-previous)), collapse=""), file = stderr()) #allow for several asterisks if an increase in one iteration is more than 2% of progress
      counter <- 0
    }
    
    F1_selfed <- rbind(create_gamete(F11, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan),create_gamete(F12, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan))
    
    if (n>0) {
      for (j in 1:n) {
        F1_selfed <- rbind(create_gamete(F1_selfed, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan),create_gamete(F1_selfed, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan))
      }
    }
    
    if (DH) {
      gamete <- create_gamete(F1_selfed, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan)
      F1_selfed <- rbind(gamete,gamete)
    }
    markersFamily <-  colSums(F1_selfed)
    # MarkersD <- markersFamily
    # MarkersD[MarkersD!=1] <- 0
    
    MarkersD <- Compute_Q(matrix(markersFamily,ncol = ncol(F1_selfed)), phi, "Genotypic")
    
    additive_values <- c(additive_values, 
                         markersFamily[QTL_numeric]%*%markereffects$YLD[QTL_numeric])
    
    genotypic_values <- c(genotypic_values, 
                          (markersFamily[QTL_numeric]%*%markereffects$YLD[QTL_numeric] + MarkersD[QTL_numeric]%*%markereffects_d$YLD[QTL_numeric]))
  }
  
  
  

  means[["true"]][["additive"]][[rep]] <- mean(additive_values)
  means[["true"]][["genotypic"]][[rep]] <- mean(genotypic_values)
  
  
  
  #progress bar
  cat("\n", file = stderr())
  
  
  additive_values_list[[rep]] <-  additive_values
  genotypic_values_list[[rep]] <-  genotypic_values
  
  
  
  
  
  

  MarkersF1 <- rbind(colSums(F11), colSums(F12))
  rownames(MarkersF1) <- c("F11", "F12")
  
  H_parentsF1 <- list(F11=F11, F12=F12)
  
  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  out$sd
  sd_hap <- c(sd_hap, out$sd)
  #Means are the same for any LD method!
  means[["MateR"]][["additive"]][[rep]] <- out$mean
  
  
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  out$sd
  sd_hap_hat <- c(sd_hap_hat, out$sd)
  #Means are the same for any LD method!
  means[["MateR_hat"]][["additive"]][[rep]] <- out$mean
  
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  out$sd
  sd_hap_hat_0.5 <- c(sd_hap_hat_0.5, out$sd)
  #Means are the same for any LD method!
  means[["MateR_hat_0.5"]][["additive"]][[rep]] <- out$mean
  
  
  
  
  
  
  #With Sigma
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  out$sd
  sd_sigma <- c(sd_sigma, out$sd)

  
  
  
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  out$sd
  sd_sigma_hat <- c(sd_sigma_hat, out$sd)

  
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  out$sd
  sd_sigma_hat_0.5 <- c(sd_sigma_hat_0.5, out$sd)

  
  
  
  
  
  
  #No LD
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind <- c(sd_ind, out$sd)

  
  
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_hat <- c(sd_ind_hat, out$sd)

  
  
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_hat_0.5 <- c(sd_ind_hat_0.5, out$sd)


  
  #With dominance

  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = markereffects_d$YLD[QTL],
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_hap_dom <- c(sd_hap_dom, out$sd)
  #Means are the same for any LD method!
  means[["MateR"]][["genotypic"]][[rep]] <- out$mean
  
  
  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = markereffects_d_hat[QTL],
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_hap_dom_hat <- c(sd_hap_dom_hat, out$sd)
  #Means are the same for any LD method!
  means[["MateR_hat"]][["genotypic"]][[rep]] <- out$mean
  
  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = markereffects_d_hat_0.5[QTL],
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_hap_dom_hat_0.5 <- c(sd_hap_dom_hat_0.5, out$sd)
  #Means are the same for any LD method!
  means[["MateR_hat_0.5"]][["genotypic"]][[rep]] <- out$mean
  

  
  
  
  
  
  
  #With Sigma
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = markereffects_d$YLD[QTL],
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = markereffects_d$YLD[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = markereffects_d$YLD[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_sigma_dom <- c(sd_sigma_dom, out$sd)

  
  
  #With Sigma
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = markereffects_d_hat[QTL],
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = markereffects_d_hat[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = markereffects_d_hat[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_sigma_dom_hat <- c(sd_sigma_dom_hat, out$sd)
  
  
  #With Sigma
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = markereffects_d_hat_0.5[QTL],
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = markereffects_d_hat_0.5[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = markereffects_d_hat_0.5[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% QTL,colnames(x) %in% QTL]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_sigma_dom_hat_0.5 <- c(sd_sigma_dom_hat_0.5, out$sd)

  
  
  
  #No LD
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = markereffects_d$YLD[QTL],
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_dom <- c(sd_ind_dom, out$sd)

  
  #No LD
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = markereffects_d_hat[QTL],
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_dom_hat <- c(sd_ind_dom_hat, out$sd)

  
  #No LD
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = markereffects_d_hat_0.5[QTL],
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_dom_hat_0.5 <- c(sd_ind_dom_hat_0.5, out$sd)
  
  
  
  ################################################################################
  #realistic scenario
  
  #additive
  
  
  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,Marker_positions],
                               phi=2,
                               markereffects=markereffects_realistic,
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,Marker_positions]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% Marker_positions,colnames(x) %in% Marker_positions]}),
                               ncores = 1)
  
  out$sd
  sd_hap_realistic <- c(sd_hap_realistic, out$sd)
  #Means are the same for any LD method!
  means[["MateR_realistic"]][["additive"]][[rep]] <- out$mean
  
  
  
  
  #With Sigma
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,Marker_positions],
                               phi=2,
                               markereffects=markereffects_realistic,
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,Marker_positions],
                                                     parametrization = "Genotypic", phi = 2, markereffects=markereffects_realistic, 
                                                     markereffects_d = NULL, 
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               
                               Sigma2=create_LD_list(Markers[,Marker_positions], 
                                                     parametrization = "Genotypic", phi = 2, 
                                                     markereffects=markereffects_realistic,
                                                     markereffects_d = NULL,
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% Marker_positions,colnames(x) %in% Marker_positions]}),
                               ncores = 1)
  
  
  out$sd
  sd_sigma_realistic <- c(sd_sigma_realistic, out$sd)
  
  
  
  #No LD
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,Marker_positions],
                               phi=2,
                               markereffects=markereffects_realistic,
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_realistic <- c(sd_ind_realistic, out$sd)
  
  
  
  
  
  
  #with dominance
  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,Marker_positions],
                               phi=2,
                               markereffects=markereffects_realistic,
                               markereffects_d = markereffects_d_realistic,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,Marker_positions]}),
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% Marker_positions,colnames(x) %in% Marker_positions]}),
                               ncores = 1)
  
  
  
  out$sd
  sd_hap_dom_realistic <- c(sd_hap_dom_realistic, out$sd)
  #Means are the same for any LD method!
  means[["MateR_realistic"]][["genotypic"]][[rep]] <- out$mean
  
  
  
  
  
  
  
  #With Sigma
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,Marker_positions],
                               phi=2,
                               markereffects=markereffects_realistic,
                               markereffects_d = markereffects_d_realistic,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,Marker_positions],
                                                     parametrization = "Genotypic", phi = 2, markereffects=markereffects_realistic, 
                                                     markereffects_d = markereffects_d_realistic, 
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               
                               Sigma2=create_LD_list(Markers[,Marker_positions], 
                                                     parametrization = "Genotypic", phi = 2, 
                                                     markereffects=markereffects_realistic,
                                                     markereffects_d = markereffects_d_realistic,
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = lapply(c_list, function(x){x[rownames(x) %in% Marker_positions,colnames(x) %in% Marker_positions]}),
                               ncores = 1)
  
  
  out$sd
  sd_sigma_dom_realistic <- c(sd_sigma_dom_realistic, out$sd)
  
  
  
  #No LD
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,Marker_positions],
                               phi=2,
                               markereffects=markereffects_realistic,
                               markereffects_d = markereffects_d_realistic,
                               n=n,
                               DHs = DH,
                               LD="Ind",
                               H_parents = NULL,
                               H_testers = NULL,
                               c_list = NULL,
                               ncores = 1)
  
  
  out$sd
  sd_ind_dom_realistic <- c(sd_ind_dom_realistic, out$sd)
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  #SimpleMating
  plan_SM <- data.frame(Parent1="F11",Parent2="F12")
  generic_Phasedgeno_SM <- rbind(H_parentsF1[[1]], H_parentsF1[[2]])
  rownames(generic_Phasedgeno_SM) <- c("F11_1", "F11_2", "F12_1", "F12_2")
  Z <- scale(MarkersF1, center = T, scale = F)
  relMat <- (Z %*% t(Z))
  relMat <- relMat/mean(diag(relMat))
  


  
  # #Simple mating does not support haplotypes for only additive!!
  # out <- getUsefA(MatePlan = plan_SM,
  #                 Markers = MarkersF1,
  #                 addEff = markereffects$YLD,
  #                 Map.In = generic_GenMap_SM,
  #                 K = relMat, #not relevant for sd computation
  #                 propSel = 0.05, #not relevant for sd computation
  #                 Type = "DH",
  #                 Generation = n)
  # 
  # sd_sigma_SM <- c(sd_sigma_SM, out[[1]]$sd)
  
  
  Criterion <- data.frame(id = c("F11", "F12"),
                          BLUPs = c(MarkersF1["F11",QTL]%*%markereffects$YLD[QTL],
                                    MarkersF1["F12",QTL]%*%markereffects$YLD[QTL]))
  means[["SM"]][["additive"]][[rep]] <- getMPV(MatePlan = plan_SM,
                                               Criterion=Criterion,
                                               K =relMat)$Y
  
  Criterion <- data.frame(id = c("F11", "F12"),
                          BLUPs = c(MarkersF1["F11",QTL]%*%markereffects_hat[QTL],
                                    MarkersF1["F12",QTL]%*%markereffects_hat[QTL]))
  means[["SM_hat"]][["additive"]][[rep]] <- getMPV(MatePlan = plan_SM,
                                               Criterion=Criterion,
                                               K =relMat)$Y
  
  Criterion <- data.frame(id = c("F11", "F12"),
                          BLUPs = c(MarkersF1["F11",QTL]%*%markereffects_hat_0.5[QTL],
                                    MarkersF1["F12",QTL]%*%markereffects_hat_0.5[QTL]))
  means[["SM_hat_0.5"]][["additive"]][[rep]] <- getMPV(MatePlan = plan_SM,
                                               Criterion=Criterion,
                                               K =relMat)$Y
  

  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = generic_Phasedgeno_SM[,QTL],
                   addEff = markereffects$YLD[QTL],
                   domEff = markereffects_d$YLD[QTL],
                   Map.In = generic_GenMap_SM[match(QTL,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "Phased")
  
  
  sd_hap_dom_SM <- c(sd_hap_dom_SM, (out[[1]]$sdA+out[[1]]$sdD))
  sd_hap_SM <- c(sd_hap_SM, out[[1]]$sdA)
  #means are the same with or without phasing!
  means[["SM"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = generic_Phasedgeno_SM[,QTL],
                   addEff = markereffects_hat[QTL],
                   domEff = markereffects_d_hat[QTL],
                   Map.In = generic_GenMap_SM[match(QTL,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "Phased")
  
  
  sd_hap_dom_SM_hat <- c(sd_hap_dom_SM_hat, (out[[1]]$sdA+out[[1]]$sdD))
  sd_hap_SM_hat <- c(sd_hap_SM_hat, out[[1]]$sdA)
  #means are the same with or without phasing!
  means[["SM_hat"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = generic_Phasedgeno_SM[,QTL],
                   addEff = markereffects_hat_0.5[QTL],
                   domEff = markereffects_d_hat_0.5[QTL],
                   Map.In = generic_GenMap_SM[match(QTL,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "Phased")
  
  
  sd_hap_dom_SM_hat_0.5 <- c(sd_hap_dom_SM_hat_0.5, (out[[1]]$sdA+out[[1]]$sdD))
  sd_hap_SM_hat_0.5 <- c(sd_hap_SM_hat_0.5, out[[1]]$sdA)
  #means are the same with or without phasing!
  means[["SM_hat_0.5"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  
  
  
  
  
  
  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = MarkersF1[,QTL],
                   addEff = markereffects$YLD[QTL],
                   domEff = markereffects_d$YLD[QTL],
                   Map.In = generic_GenMap_SM[match(QTL,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "NonPhased")
  
  
  sd_sigma_dom_SM <- c(sd_sigma_dom_SM, (out[[1]]$sdA+out[[1]]$sdD))
  sd_sigma_SM <- c(sd_sigma_SM, out[[1]]$sdA)
 # means[["sigma_SM"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = MarkersF1[,QTL],
                   addEff = markereffects_hat[QTL],
                   domEff = markereffects_d_hat[QTL],
                   Map.In = generic_GenMap_SM[match(QTL,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "NonPhased")
  
  
  sd_sigma_dom_SM_hat <- c(sd_sigma_dom_SM_hat, (out[[1]]$sdA+out[[1]]$sdD))
  sd_sigma_SM_hat <- c(sd_sigma_SM_hat, out[[1]]$sdA)
 # means[["sigma_SM_hat"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = MarkersF1[,QTL],
                   addEff = markereffects_hat_0.5[QTL],
                   domEff = markereffects_d_hat_0.5[QTL],
                   Map.In = generic_GenMap_SM[match(QTL,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "NonPhased")
  
  
  sd_sigma_dom_SM_hat_0.5 <- c(sd_sigma_dom_SM_hat_0.5, (out[[1]]$sdA+out[[1]]$sdD))
  sd_sigma_SM_hat_0.5 <- c(sd_sigma_SM_hat_0.5, out[[1]]$sdA)
#  means[["sigma_SM_hat_0.5"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  
  
  
  
  ################################################################################
  #realistic scenario
  
  Criterion <- data.frame(id = c("F11", "F12"),
                          BLUPs = c(MarkersF1["F11",][Marker_positions]%*%markereffects_realistic,
                                    MarkersF1["F12",][Marker_positions]%*%markereffects_realistic))
  means[["SM_realistic"]][["additive"]][[rep]] <- getMPV(MatePlan = plan_SM,
                                                         Criterion=Criterion,
                                                         K =relMat)$Y
  
  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = generic_Phasedgeno_SM[,Marker_positions],
                   addEff = markereffects_realistic,
                   domEff = markereffects_d_realistic,
                   Map.In = generic_GenMap_SM[match(Marker_positions,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "Phased")
  
  
  sd_hap_dom_SM_realistic <- c(sd_hap_dom_SM_realistic, (out[[1]]$sdA+out[[1]]$sdD))
  sd_hap_SM_realistic <- c(sd_hap_SM_realistic, out[[1]]$sdA)
  #means are the same with or without phasing!
  means[["SM_realistic"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  
  
  out <- getUsefAD(MatePlan = plan_SM,
                   Markers = MarkersF1[,Marker_positions],
                   addEff = markereffects_realistic,
                   domEff = markereffects_d_realistic,
                   Map.In = generic_GenMap_SM[match(Marker_positions,generic_GenMap_SM$mkr),],
                   K = relMat, #not relevant for sd computation
                   propSel = 0.05, #not relevant for sd computation
                   Method = "NonPhased")
  
  
  sd_sigma_dom_SM_realistic <- c(sd_sigma_dom_SM_realistic, (out[[1]]$sdA+out[[1]]$sdD))
  sd_sigma_SM_realistic <- c(sd_sigma_SM_realistic, out[[1]]$sdA)
  # means[["sigma_SM"]][["genotypic"]][[rep]] <- out[[1]]$Total.gv
  
  
  
  
  
  
  
  
  
  #wolfe 2021 (predCrossVar) implementation
  
  
  c_mat <- matrix(0, nrow = length(markereffects$YLD), ncol = length(markereffects$YLD))
  counter <- 1
  for (i in 1:length(c_list)) {
    tmp <- c_list[[i]]
    c_mat[counter:(counter+nrow(tmp)-1),counter:(counter+nrow(tmp)-1)] <- tmp
    counter <- counter + nrow(tmp)
  }
  recombFreqMat <- 1-2*c_mat
  rownames(recombFreqMat) <- colnames(recombFreqMat) <- colnames(Markers)
  
  haploMat <- rbind(H_parentsF1[[1]], H_parentsF1[[2]])
  rownames(haploMat) <- c("F11_HapA","F11_HapB","F12_HapA","F12_HapB")
  
  
  addEffectsT1 <- matrix(markereffects$YLD[QTL], ncol=1)
  rownames(addEffectsT1) <- colnames(haploMat[,QTL])
  domEffectsT1<-matrix(markereffects_d$YLD[QTL], ncol=1)
  rownames(domEffectsT1) <- colnames(haploMat[,QTL])
  
  ped2predict<-crosses2predict(parents = c("F11","F12"))
  #ped2predict
  #remove selfing F11 or selfing F12, leave only F11xF12
  ped2predict <- ped2predict[2,]
  
  
  out <- predCrossVarA(sireID = "F11", damID = "F12", 
                addEffects = addEffectsT1,
                haploMat = haploMat[,QTL],
                recombFreqMat = recombFreqMat[QTL,QTL])

  sd_wolfe <- c(sd_wolfe, sqrt(c(out$varA)))
  means[["wolfe"]][["additive"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                          damID="F12",
                                                          doseMat=MarkersF1[,QTL],
                                                          postMeanAddEffects=addEffectsT1,
                                                          postMeanDomEffects=NULL,
                                                          q_pop = NULL,
                                                          p_pop = NULL)
  
  

  
  
  out <- runCrossVarPredsAD(ped = ped2predict,
                            addEffects = addEffectsT1,
                            domEffects = domEffectsT1,
                            haploMat = haploMat[,QTL],
                            recombFreqMat = recombFreqMat[QTL,QTL])

  sd_wolfe_dom <- c(sd_wolfe_dom, sqrt(c(out$varA+out$varD)))
  means[["wolfe"]][["genotypic"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                      damID="F12",
                                                      doseMat=MarkersF1[,QTL],
                                                      postMeanAddEffects=addEffectsT1,
                                                      postMeanDomEffects=domEffectsT1,
                                                      q_pop = (1-colMeans(MarkersF1)/phi),
                                                      p_pop = colMeans(MarkersF1)/phi)
  
  
  
  
  
  
  
  addEffectsT1 <- matrix(markereffects_hat[QTL], ncol=1)
  rownames(addEffectsT1) <- colnames(haploMat[,QTL])
  domEffectsT1<-matrix(markereffects_d_hat[QTL], ncol=1)
  rownames(domEffectsT1) <- colnames(haploMat[,QTL])
  
  
  out <- predCrossVarA(sireID = "F11", damID = "F12", 
                       addEffects = addEffectsT1,
                       haploMat = haploMat[,QTL],
                       recombFreqMat = recombFreqMat[QTL,QTL])
  
  sd_wolfe_hat <- c(sd_wolfe_hat, sqrt(c(out$varA)))
  means[["wolfe_hat"]][["additive"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                        damID="F12",
                                                        doseMat=MarkersF1[,QTL],
                                                        postMeanAddEffects=addEffectsT1,
                                                        postMeanDomEffects=NULL,
                                                        q_pop = NULL,
                                                        p_pop = NULL)

  
  out <- runCrossVarPredsAD(ped = ped2predict,
                            addEffects = addEffectsT1,
                            domEffects = domEffectsT1,
                            haploMat = haploMat[,QTL],
                            recombFreqMat = recombFreqMat[QTL,QTL])
  
  sd_wolfe_dom_hat <- c(sd_wolfe_dom_hat, sqrt(c(out$varA+out$varD)))
  means[["wolfe_hat"]][["genotypic"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                         damID="F12",
                                                         doseMat=MarkersF1[,QTL],
                                                         postMeanAddEffects=addEffectsT1,
                                                         postMeanDomEffects=domEffectsT1,
                                                         q_pop = (1-colMeans(MarkersF1[,QTL])/phi),
                                                         p_pop = colMeans(MarkersF1[,QTL])/phi)
  
  
  
  
  addEffectsT1 <- matrix(markereffects_hat_0.5[QTL], ncol=1)
  rownames(addEffectsT1) <- colnames(haploMat[,QTL])
  domEffectsT1<-matrix(markereffects_d_hat_0.5[QTL], ncol=1)
  rownames(domEffectsT1) <- colnames(haploMat[,QTL])
  

  
  
  out <- predCrossVarA(sireID = "F11", damID = "F12", 
                       addEffects = addEffectsT1,
                       haploMat = haploMat[,QTL],
                       recombFreqMat = recombFreqMat[QTL,QTL])
  
  sd_wolfe_hat_0.5 <- c(sd_wolfe_hat_0.5, sqrt(c(out$varA)))
  means[["wolfe_hat_0.5"]][["additive"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                        damID="F12",
                                                        doseMat=MarkersF1[,QTL],
                                                        postMeanAddEffects=addEffectsT1,
                                                        postMeanDomEffects=NULL,
                                                        q_pop = NULL,
                                                        p_pop = NULL)
  
  
  out <- runCrossVarPredsAD(ped = ped2predict,
                            addEffects = addEffectsT1,
                            domEffects = domEffectsT1,
                            haploMat = haploMat[,QTL],
                            recombFreqMat = recombFreqMat[QTL,QTL])
  
  sd_wolfe_dom_hat_0.5 <- c(sd_wolfe_dom_hat_0.5, sqrt(c(out$varA+out$varD)))
  means[["wolfe_hat_0.5"]][["genotypic"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                         damID="F12",
                                                         doseMat=MarkersF1[,QTL],
                                                         postMeanAddEffects=addEffectsT1,
                                                         postMeanDomEffects=domEffectsT1,
                                                         q_pop = (1-colMeans(MarkersF1[,QTL])/phi),
                                                         p_pop = colMeans(MarkersF1[,QTL])/phi)
  
  
  
  ################################################################################
  #realistic scenario
  
  
  addEffectsT1 <- matrix(markereffects_realistic, ncol=1)
  rownames(addEffectsT1) <- colnames(haploMat[,Marker_positions])
  domEffectsT1<-matrix(markereffects_d_realistic, ncol=1)
  rownames(domEffectsT1) <- colnames(haploMat[,Marker_positions])
  
  
  
  
  out <- predCrossVarA(sireID = "F11", damID = "F12", 
                       addEffects = addEffectsT1,
                       haploMat = haploMat[,Marker_positions],
                       recombFreqMat = recombFreqMat[Marker_positions,Marker_positions])
  
  sd_wolfe_realistic <- c(sd_wolfe_realistic, sqrt(c(out$varA)))
  means[["wolfe_realistic"]][["additive"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                               damID="F12",
                                                               doseMat=MarkersF1[,Marker_positions],
                                                               postMeanAddEffects=addEffectsT1,
                                                               postMeanDomEffects=NULL,
                                                               q_pop = NULL,
                                                               p_pop = NULL)
  
  
  out <- runCrossVarPredsAD(ped = ped2predict,
                            addEffects = addEffectsT1,
                            domEffects = domEffectsT1,
                            haploMat = haploMat[,Marker_positions],
                            recombFreqMat = recombFreqMat[Marker_positions,Marker_positions])
  
  sd_wolfe_dom_realistic <- c(sd_wolfe_dom_realistic, sqrt(c(out$varA+out$varD)))
  means[["wolfe_realistic"]][["genotypic"]][[rep]] <- Wolfe_GVs(sireID="F11",
                                                                damID="F12",
                                                                doseMat=MarkersF1[,Marker_positions],
                                                                postMeanAddEffects=addEffectsT1,
                                                                postMeanDomEffects=domEffectsT1,
                                                                q_pop = (1-colMeans(MarkersF1[,Marker_positions])/phi),
                                                                p_pop = colMeans(MarkersF1[,Marker_positions])/phi)
  
  
  
  
  
  
  
  
  ########
  #GenomicMating
  MarkersGM <- MarkersF1-1
  Z <- matrix(0.5, nrow = 1, ncol = 2)
  
  out <- getstatsM1R(MARKERS = MarkersGM[,QTL],markereffects = markereffects$YLD[QTL],P = Z)
  sd_GM <- c(sd_GM, sqrt(out[3]))
  means[["GM"]][["additive"]][[rep]] <- out[2]
  
  


  out <- getstatsM1R(MARKERS = MarkersGM[,QTL],markereffects = markereffects_hat[QTL] ,P = Z)
  sd_GM_hat <- c(sd_GM_hat, sqrt(out[3]))
  means[["GM_hat"]][["additive"]][[rep]] <- out[2]
  


  out <- getstatsM1R(MARKERS = MarkersGM[,QTL],markereffects = markereffects_hat_0.5[QTL],P = Z)
  sd_GM_hat_0.5 <- c(sd_GM_hat_0.5, sqrt(out[3]))
  means[["GM_hat_0.5"]][["additive"]][[rep]] <- out[2]



  ################################################################################
  #realistic scenario
  
  
  out <- getstatsM1R(MARKERS = MarkersGM[,Marker_positions],
                     markereffects = markereffects_realistic,
                     P = Z)
  sd_GM_realistic <- c(sd_GM_realistic, sqrt(out[3]))
  means[["GM_realistic"]][["additive"]][[rep]] <- out[2]
  
  
  
  
  
  
}

results_clonal <- list(additive_values_list=additive_values_list,
                   genotypic_values_list=genotypic_values_list,
                   parents=parents,
                   means=means,
                   
                   sd_hap=sd_hap,
                   sd_sigma=sd_sigma,
                   sd_ind=sd_ind,
                   sd_hap_dom=sd_hap_dom,
                   sd_sigma_dom=sd_sigma_dom,
                   sd_ind_dom=sd_ind_dom,
                   sd_wolfe = sd_wolfe,
                   sd_wolfe_dom=sd_wolfe_dom,
                   sd_hap_SM=sd_hap_SM,
                   sd_sigma_SM=sd_sigma_SM,
                   sd_hap_dom_SM=sd_hap_dom_SM,
                   sd_sigma_dom_SM=sd_sigma_dom_SM
                   
                   ,
                   sd_hap_hat=sd_hap_hat,
                   sd_sigma_hat=sd_sigma_hat,
                   sd_ind_hat=sd_ind_hat,
                   sd_hap_dom_hat=sd_hap_dom_hat,
                   sd_sigma_dom_hat=sd_sigma_dom_hat,
                   sd_ind_dom_hat=sd_ind_dom_hat,
                   sd_wolfe_hat = sd_wolfe_hat,
                   sd_wolfe_dom_hat=sd_wolfe_dom_hat,
                   sd_hap_SM_hat=sd_hap_SM_hat,
                   sd_sigma_SM_hat=sd_sigma_SM_hat,
                   sd_hap_dom_SM_hat=sd_hap_dom_SM_hat,
                   sd_sigma_dom_SM_hat=sd_sigma_dom_SM_hat
                   
                   ,
                   sd_hap_hat_0.5=sd_hap_hat_0.5,
                   sd_sigma_hat_0.5=sd_sigma_hat_0.5,
                   sd_ind_hat_0.5=sd_ind_hat_0.5,
                   sd_hap_dom_hat_0.5=sd_hap_dom_hat_0.5,
                   sd_sigma_dom_hat_0.5=sd_sigma_dom_hat_0.5,
                   sd_ind_dom_hat_0.5=sd_ind_dom_hat_0.5,
                   sd_wolfe_hat_0.5 = sd_wolfe_hat_0.5,
                   sd_wolfe_dom_hat_0.5=sd_wolfe_dom_hat_0.5,
                   sd_hap_SM_hat_0.5=sd_hap_SM_hat_0.5,
                   sd_sigma_SM_hat_0.5=sd_sigma_SM_hat_0.5,
                   sd_hap_dom_SM_hat_0.5=sd_hap_dom_SM_hat_0.5,
                   sd_sigma_dom_SM_hat_0.5=sd_sigma_dom_SM_hat_0.5,
                   
                   sd_hap_realistic=sd_hap_realistic,
                   sd_sigma_realistic=sd_sigma_realistic,
                   sd_ind_realistic=sd_ind_realistic,
                   sd_hap_dom_realistic=sd_hap_dom_realistic,
                   sd_sigma_dom_realistic=sd_sigma_dom_realistic,
                   sd_ind_dom_realistic=sd_ind_dom_realistic,
                   
                   #wolfe 2021
                   sd_wolfe_realistic=sd_wolfe_realistic,
                   sd_wolfe_dom_realistic=sd_wolfe_dom_realistic,
                   
                   #Simple Mating
                   sd_hap_SM_realistic=sd_hap_SM_realistic,
                   sd_sigma_SM_realistic=sd_sigma_SM_realistic,
                   sd_hap_dom_SM_realistic=sd_hap_dom_SM_realistic,
                   sd_sigma_dom_SM_realistic=sd_sigma_dom_SM_realistic,

                   sd_GM=sd_GM,
                   sd_GM_hat=sd_GM_hat,
                   sd_GM_hat_0.5=sd_GM_hat_0.5,
                   sd_GM_realistic=sd_GM_realistic)


save(results_clonal, file = "./Results_clonal.RData")



sd_additive <- c()
sd_genotypic <- c()
for (i in 1:nrep) {
  sd_additive <- c(sd_additive, sd(additive_values_list[[i]]))
  sd_genotypic <- c(sd_genotypic, sd(genotypic_values_list[[i]]))
}

sd_additive
sd_genotypic



#correlations and RMSE:

################################################################################
#True Marker Effects
################################################################################

#########################
#only additive
#########################

#correlations
#MateR
cor(sd_additive, unlist(sd_hap))
cor(sd_additive, unlist(sd_sigma))
cor(sd_additive, unlist(sd_ind))

#SimpleMating
cor(sd_additive, sd_hap_SM)
cor(sd_additive, sd_sigma_SM)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe[is.na(sd_wolfe)] <- 0 
cor(sd_additive, sd_wolfe)


#GenomicMating
cor(sd_additive, sd_GM)



#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind))^2 ))/mean(sd_additive)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_SM)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_SM)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe)^2 ))/mean(sd_genotypic)*100


#GenomicMating
sqrt(mean( (sd_genotypic-sd_GM)^2 ))/mean(sd_genotypic)*100


#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom))
cor(sd_genotypic, unlist(sd_sigma_dom))
cor(sd_genotypic, unlist(sd_ind_dom))

#SimpleMating
cor(sd_additive, sd_hap_dom_SM)
cor(sd_additive, sd_sigma_dom_SM)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_dom[is.na(sd_wolfe_dom)] <- 0 
cor(sd_additive, sd_wolfe_dom)



#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom))^2 ))/mean(sd_genotypic)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_dom_SM)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_dom_SM)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_dom)^2 ))/mean(sd_genotypic)*100







################################################################################
#Estimated Marker Effects, GS accuracy = 1
################################################################################

#########################
#only additive
#########################

#correlations
#MateR
cor(sd_additive, unlist(sd_hap_hat))
cor(sd_additive, unlist(sd_sigma_hat))
cor(sd_additive, unlist(sd_ind_hat))

#SimpleMating
cor(sd_additive, sd_hap_SM_hat)
cor(sd_additive, sd_sigma_SM_hat)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_hat[is.na(sd_wolfe_hat)] <- 0 
cor(sd_additive, sd_wolfe_hat)

#GenomicMating
cor(sd_additive, sd_GM_hat)




#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap_hat))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma_hat))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind_hat))^2 ))/mean(sd_additive)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_SM_hat)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_SM_hat)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_hat)^2 ))/mean(sd_genotypic)*100

#GenomicMating
sqrt(mean( (sd_genotypic-sd_GM_hat)^2 ))/mean(sd_genotypic)*100

#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom_hat))
cor(sd_genotypic, unlist(sd_sigma_dom_hat))
cor(sd_genotypic, unlist(sd_ind_dom_hat))

#SimpleMating
cor(sd_additive, sd_hap_dom_SM_hat)
cor(sd_additive, sd_sigma_dom_SM_hat)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_dom_hat[is.na(sd_wolfe_dom_hat)] <- 0 
cor(sd_additive, sd_wolfe_dom_hat)



#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom_hat))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom_hat))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom_hat))^2 ))/mean(sd_genotypic)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_dom_SM_hat)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_dom_SM_hat)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_dom_hat)^2 ))/mean(sd_genotypic)*100











################################################################################
#Estimated Marker Effects, GS accuracy = sqrt(0.5)
################################################################################

#########################
#only additive
#########################

#correlations
#MateR
cor(sd_additive, unlist(sd_hap_hat_0.5))
cor(sd_additive, unlist(sd_sigma_hat_0.5))
cor(sd_additive, unlist(sd_ind_hat_0.5))

#SimpleMating
cor(sd_additive, sd_hap_SM_hat_0.5)
cor(sd_additive, sd_sigma_SM_hat_0.5)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_hat_0.5[is.na(sd_wolfe_hat_0.5)] <- 0 
cor(sd_additive, sd_wolfe_hat_0.5)

#GenomicMating
cor(sd_additive, sd_GM_hat_0.5)




#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap_hat_0.5))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma_hat_0.5))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind_hat_0.5))^2 ))/mean(sd_additive)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_SM_hat_0.5)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_SM_hat_0.5)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_hat_0.5)^2 ))/mean(sd_genotypic)*100

#GenomicMating
sqrt(mean( (sd_genotypic-sd_GM_hat_0.5)^2 ))/mean(sd_genotypic)*100

#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom_hat_0.5))
cor(sd_genotypic, unlist(sd_sigma_dom_hat_0.5))
cor(sd_genotypic, unlist(sd_ind_dom_hat_0.5))

#SimpleMating
cor(sd_additive, sd_hap_dom_SM_hat_0.5)
cor(sd_additive, sd_sigma_dom_SM_hat_0.5)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_dom_hat_0.5[is.na(sd_wolfe_dom_hat_0.5)] <- 0 
cor(sd_additive, sd_wolfe_dom_hat_0.5)



#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom_hat_0.5))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom_hat_0.5))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom_hat_0.5))^2 ))/mean(sd_genotypic)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_dom_SM_hat_0.5)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_dom_SM_hat_0.5)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_dom_hat_0.5)^2 ))/mean(sd_genotypic)*100





################################################################################
#Markers != QTLs
#Estimated Marker Effects, GS accuracy = sqrt(0.5)
################################################################################

#########################
#only additive
#########################

#correlations
#MateR
cor(sd_additive, unlist(sd_hap_realistic))
cor(sd_additive, unlist(sd_sigma_realistic))
cor(sd_additive, unlist(sd_ind_realistic))

#SimpleMating
cor(sd_additive, sd_hap_SM_realistic)
cor(sd_additive, sd_sigma_SM_realistic)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_realistic[is.na(sd_wolfe_realistic)] <- 0 
cor(sd_additive, sd_wolfe_realistic)

#GenomicMating
cor(sd_additive, sd_GM_realistic)




#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap_realistic))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma_realistic))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind_realistic))^2 ))/mean(sd_additive)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_SM_realistic)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_SM_realistic)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_realistic)^2 ))/mean(sd_genotypic)*100

#GenomicMating
sqrt(mean( (sd_genotypic-sd_GM_realistic)^2 ))/mean(sd_genotypic)*100

#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom_realistic))
cor(sd_genotypic, unlist(sd_sigma_dom_realistic))
cor(sd_genotypic, unlist(sd_ind_dom_realistic))

#SimpleMating
cor(sd_additive, sd_hap_dom_SM_realistic)
cor(sd_additive, sd_sigma_dom_SM_realistic)

#Wolfe
#Negative variances cause NaN sd. Set that to zero
sd_wolfe_dom_realistic[is.na(sd_wolfe_dom_realistic)] <- 0 
cor(sd_additive, sd_wolfe_dom_realistic)



#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom_realistic))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom_realistic))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom_realistic))^2 ))/mean(sd_genotypic)*100

#SimpleMating
sqrt(mean( (sd_genotypic-sd_hap_dom_SM_realistic)^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-sd_sigma_dom_SM_realistic)^2 ))/mean(sd_genotypic)*100

#wolfe 
sqrt(mean( (sd_genotypic-sd_wolfe_dom_realistic)^2 ))/mean(sd_genotypic)*100








###############################################################################
#means
###############################################################################

#additive

#MateR
cor(unlist(means$true$additive), unlist(means$MateR$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$MateR$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$MateR_hat$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$MateR_hat$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$MateR_hat_0.5$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$MateR_hat_0.5$additive))^2 ))/mean(unlist(means$true$additive))*100


cor(unlist(means$true$additive), unlist(means$MateR_realistic$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$MateR_realistic$additive))^2 ))/mean(unlist(means$true$additive))*100




#SimpleMating
cor(unlist(means$true$additive), unlist(means$SM$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$SM$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$SM_hat$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$SM_hat$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$SM_hat_0.5$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$SM_hat_0.5$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$SM_realistic$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$SM_realistic$additive))^2 ))/mean(unlist(means$true$additive))*100



#Wolfe
cor(unlist(means$true$additive), unlist(means$wolfe$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$wolfe$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$wolfe_hat$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$wolfe_hat$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$wolfe_hat_0.5$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$wolfe_hat_0.5$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$wolfe_realistic$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$wolfe_realistic$additive))^2 ))/mean(unlist(means$true$additive))*100




#GenomicMating
cor(unlist(means$true$additive), unlist(means$GM$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$GM$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$GM_hat$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$GM_hat$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$GM_hat_0.5$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$GM_hat_0.5$additive))^2 ))/mean(unlist(means$true$additive))*100

cor(unlist(means$true$additive), unlist(means$GM_realistic$additive))
sqrt(mean( (unlist(means$true$additive)-unlist(means$GM_realistic$additive))^2 ))/mean(unlist(means$true$additive))*100




#Genotypic

#MateR
cor(unlist(means$true$genotypic), unlist(means$MateR$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$MateR$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$MateR_hat$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$MateR_hat$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$MateR_hat_0.5$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$MateR_hat_0.5$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$MateR_realistic$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$MateR_realistic$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100




#SimpleMating
cor(unlist(means$true$genotypic), unlist(means$SM$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$SM$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$SM_hat$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$SM_hat$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$SM_hat_0.5$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$SM_hat_0.5$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$SM_realistic$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$SM_realistic$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100




#Wolfe
cor(unlist(means$true$genotypic), unlist(means$wolfe$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$wolfe$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$wolfe_hat$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$wolfe_hat$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$wolfe_hat_0.5$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$wolfe_hat_0.5$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

cor(unlist(means$true$genotypic), unlist(means$wolfe_realistic$genotypic))
sqrt(mean( (unlist(means$true$genotypic)-unlist(means$wolfe_realistic$genotypic))^2 ))/mean(unlist(means$true$genotypic))*100

