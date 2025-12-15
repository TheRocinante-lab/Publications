set.seed(1234)
n <- 0
n2 <- 10000
nrep <- 500
DH=T
phi = 2
nQTL = 250


library(MateR)
library(SimpleMating)
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




#Homozygous parents, DHs

#Initialize vectors for results of within-family sd.
#In hindsight, I should have just put them in a list instead of initializing
#so many vectors. I did that for the family mean. 

#naming convention:
#the addition of "hap" indicates that full haplotype data was used to compute LD
#the addition of "sigma" indicates LD was approximated from the parental population
#the addition of "ind" indicates LD was disregarded
#the addition of "dom" indicates dominance effects were considered


additive_values_list <-  list()
genotypic_values_list <-  list()
parents <- list()
testers <- list()

#MateR
sd_hap <- c()
sd_sigma <- c()
sd_ind <- c()
sd_hap_dom <- c()
sd_sigma_dom <- c()
sd_ind_dom <- c()


#MateR
sd_hap_hat <- c()
sd_sigma_hat <- c()
sd_ind_hat <- c()
sd_hap_dom_hat <- c()
sd_sigma_dom_hat <- c()
sd_ind_dom_hat <- c()

#MateR
sd_hap_hat_0.5 <- c()
sd_sigma_hat_0.5 <- c()
sd_ind_hat_0.5 <- c()
sd_hap_dom_hat_0.5 <- c()
sd_sigma_dom_hat_0.5 <- c()
sd_ind_dom_hat_0.5 <- c()

#MateR
sd_hap_realistic <- c()
sd_sigma_realistic <- c()
sd_ind_realistic <- c()
sd_hap_dom_realistic <- c()
sd_sigma_dom_realistic <- c()
sd_ind_dom_realistic <- c()

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
  
  
  #Use directly the parental lines here, not the F1 from them!
  parents1 = "tmp"
  parents2 = "tmp" 
  while (parents1 == parents2) {
    parents1 <- sample(rownames(Markers),1)
    F11 <- H_parents[[parents1]]
    parents2 <- sample(rownames(Markers),1)
    F12 <- H_parents[[parents2]]
  }
  
  parents[[rep]] <- c(parents1, parents2)
  
  tester =parents1
  while (tester %in% c(parents1,parents2)) {
    tester <- sample(rownames(Markers),1)
    Tester <- H_parents[[tester]]
  }
  
  testers[[rep]] <- tester
  
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
    
    #Cross with tester
    gamete <- create_gamete(F1_selfed, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan)
    gameteT <- create_gamete(Tester, Chromosomes, ChromosomesCentiMorgan = ChromosomesCentiMorgan)
    F1_selfed <- rbind(gamete,gameteT)
    
    markersFamily <-  colSums(F1_selfed)
    
    MarkersD <- Compute_Q(matrix(markersFamily,ncol = ncol(F1_selfed)), phi, "Genotypic")
    
    additive_values <- c(additive_values, 
                         markersFamily[QTL_numeric]%*%markereffects$YLD[QTL_numeric])
    
    genotypic_values <- c(genotypic_values, 
                          (markersFamily[QTL_numeric]%*%markereffects$YLD[QTL_numeric] + MarkersD[QTL_numeric]%*%markereffects_d$YLD[QTL_numeric]))
    
    
  }
  
  #progress bar
  cat("\n", file = stderr())
  
  
  additive_values_list[[rep]] <-  additive_values
  genotypic_values_list[[rep]] <-  genotypic_values
  
  
  
  means[["true"]][["additive"]][[rep]] <- mean(additive_values)
  means[["true"]][["genotypic"]][[rep]] <- mean(genotypic_values)
  
  
  
  MarkersF1 <- rbind(colSums(F11), colSums(F12))
  rownames(MarkersF1) <- c("F11", "F12")
  
  MarkersT <- matrix(colSums(Tester), nrow = 1)
  rownames(MarkersT) <- c("Tester")
  colnames(MarkersT) <- colnames(Markers)
  
  H_parentsF1 <- list(F11=F11, F12=F12)
  H_tester <- list(Tester = Tester)
  
  
  #With Haplotypes
  P <- make_P_matrix("F11", "F12")
  out <- calculateFamilyValues(P=P,
                               parametrization = "Genotypic",
                               Username=username, 
                               Password=password,
                               Markers=MarkersF1[,QTL],
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = lapply(H_tester, function(x){x[,QTL]}),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = lapply(H_tester, function(x){x[,QTL]}),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Full",
                               H_parents = lapply(H_parentsF1, function(x){x[,QTL]}),
                               H_testers = lapply(H_tester, function(x){x[,QTL]}),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               SigmaT=create_LD_list(rbind(Markers[,QTL], MarkersT[,QTL, drop = F]), parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               SigmaT=create_LD_list(rbind(Markers[,QTL], MarkersT[,QTL, drop = F]), parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               SigmaT=create_LD_list(rbind(Markers[,QTL], MarkersT[,QTL, drop = F]), parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = NULL, Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects$YLD[QTL],
                               markereffects_d = markereffects_d$YLD[QTL],
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = markereffects_d$YLD[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = markereffects_d$YLD[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               SigmaT=create_LD_list(rbind(Markers[,QTL], MarkersT[,QTL, drop = F]), parametrization = "Genotypic", phi = 2, markereffects=markereffects$YLD[QTL], markereffects_d = markereffects_d$YLD[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects_hat[QTL],
                               markereffects_d = markereffects_d_hat[QTL],
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = markereffects_d_hat[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = markereffects_d_hat[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               SigmaT=create_LD_list(rbind(Markers[,QTL], MarkersT[,QTL, drop = F]), parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat[QTL], markereffects_d = markereffects_d_hat[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
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
                               Markers_T = MarkersT[,QTL, drop = F],
                               phi=2,
                               markereffects=markereffects_hat_0.5[QTL],
                               markereffects_d = markereffects_d_hat_0.5[QTL],
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = markereffects_d_hat_0.5[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               Sigma2=create_LD_list(Markers[,QTL], parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = markereffects_d_hat_0.5[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
                               SigmaT=create_LD_list(rbind(Markers[,QTL], MarkersT[,QTL, drop = F]), parametrization = "Genotypic", phi = 2, markereffects=markereffects_hat_0.5[QTL], markereffects_d = markereffects_d_hat_0.5[QTL], Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% QTL]})),
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,QTL, drop = F],
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
                               Markers_T = MarkersT[,Marker_positions, drop=F],
                               phi=phi,
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
                               Markers_T = MarkersT[,Marker_positions, drop=F],
                               phi=phi,
                               markereffects=markereffects_realistic,
                               markereffects_d = NULL,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,Marker_positions],
                                                     parametrization = "Genotypic", phi = phi, markereffects=markereffects_realistic, 
                                                     markereffects_d = NULL, 
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               
                               Sigma2=create_LD_list(Markers[,Marker_positions], 
                                                     parametrization = "Genotypic", phi = phi, 
                                                     markereffects=markereffects_realistic,
                                                     markereffects_d = NULL,
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               SigmaT=create_LD_list(rbind(Markers[,Marker_positions], MarkersT[,Marker_positions]), 
                                                     parametrization = "Genotypic", phi = phi,
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
                               Markers_T = MarkersT[,Marker_positions, drop=F],
                               phi=phi,
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
                               Markers_T = MarkersT[,Marker_positions, drop=F],
                               phi=phi,
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
                               Markers_T = MarkersT[,Marker_positions, drop=F],
                               phi=phi,
                               markereffects=markereffects_realistic,
                               markereffects_d = markereffects_d_realistic,
                               n=n,
                               DHs = DH,
                               LD="Approx",
                               Sigma1=create_LD_list(Markers[,Marker_positions],
                                                     parametrization = "Genotypic", phi = phi, markereffects=markereffects_realistic, 
                                                     markereffects_d = markereffects_d_realistic, 
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               
                               Sigma2=create_LD_list(Markers[,Marker_positions], 
                                                     parametrization = "Genotypic", phi = phi, 
                                                     markereffects=markereffects_realistic,
                                                     markereffects_d = markereffects_d_realistic,
                                                     Chromosomes = lapply(Chromosomes_names, function(x){x[x %in% Marker_positions]})),
                               SigmaT=create_LD_list(rbind(Markers[,Marker_positions], MarkersT[,Marker_positions]), 
                                                     parametrization = "Genotypic", phi = phi,
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
                               Markers_T = MarkersT[,Marker_positions, drop=F],
                               phi=phi,
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
  
  
  
  
  
  
  
  
  
  
}

results_DHs_Tester <- list(additive_values_list=additive_values_list,
                            genotypic_values_list=genotypic_values_list,

                           parents=parents,
                           means=means,
                           testers=testers,
                           
                            sd_hap=sd_hap,
                            sd_sigma=sd_sigma,
                            sd_ind=sd_ind,
                            sd_hap_dom=sd_hap_dom,
                            sd_sigma_dom=sd_sigma_dom,
                            sd_ind_dom=sd_ind_dom,
                           
                           sd_hap_hat=sd_hap_hat,
                           sd_sigma_hat=sd_sigma_hat,
                           sd_ind_hat=sd_ind_hat,
                           sd_hap_dom_hat=sd_hap_dom_hat,
                           sd_sigma_dom_hat=sd_sigma_dom_hat,
                           sd_ind_dom_hat=sd_ind_dom_hat,
                           
                           sd_hap_hat_0.5=sd_hap_hat_0.5,
                           sd_sigma_hat_0.5=sd_sigma_hat_0.5,
                           sd_ind_hat_0.5=sd_ind_hat_0.5,
                           sd_hap_dom_hat_0.5=sd_hap_dom_hat_0.5,
                           sd_sigma_dom_hat_0.5=sd_sigma_dom_hat_0.5,
                           sd_ind_dom_hat_0.5=sd_ind_dom_hat_0.5,
                           
                           sd_hap_realistic=sd_hap_realistic,
                           sd_sigma_realistic=sd_sigma_realistic,
                           sd_ind_realistic=sd_ind_realistic,
                           sd_hap_dom_realistic=sd_hap_dom_realistic,
                           sd_sigma_dom_realistic=sd_sigma_dom_realistic,
                           sd_ind_dom_realistic=sd_ind_dom_realistic
)


save(results_DHs_Tester, file = paste0("./ResultsDH_Tester_n_", n,".RData"))



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



#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind))^2 ))/mean(sd_additive)*100

#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom))
cor(sd_genotypic, unlist(sd_sigma_dom))
cor(sd_genotypic, unlist(sd_ind_dom))


#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom))^2 ))/mean(sd_genotypic)*100




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


#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap_hat))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma_hat))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind_hat))^2 ))/mean(sd_additive)*100


#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom_hat))
cor(sd_genotypic, unlist(sd_sigma_dom_hat))
cor(sd_genotypic, unlist(sd_ind_dom_hat))


#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom_hat))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom_hat))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom_hat))^2 ))/mean(sd_genotypic)*100



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


#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap_hat_0.5))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma_hat_0.5))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind_hat_0.5))^2 ))/mean(sd_additive)*100



#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom_hat_0.5))
cor(sd_genotypic, unlist(sd_sigma_dom_hat_0.5))
cor(sd_genotypic, unlist(sd_ind_dom_hat_0.5))


#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom_hat_0.5))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom_hat_0.5))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom_hat_0.5))^2 ))/mean(sd_genotypic)*100







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


#normalized RMSE, as percentage
sqrt(mean( (sd_additive-unlist(sd_hap_realistic))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_sigma_realistic))^2 ))/mean(sd_additive)*100
sqrt(mean( (sd_additive-unlist(sd_ind_realistic))^2 ))/mean(sd_additive)*100


#########################
#with dominance
#########################

#correlations
#MateR
cor(sd_genotypic , unlist(sd_hap_dom_realistic))
cor(sd_genotypic, unlist(sd_sigma_dom_realistic))
cor(sd_genotypic, unlist(sd_ind_dom_realistic))



#normalized RMSE, as percentage
sqrt(mean( (sd_genotypic-unlist(sd_hap_dom_realistic))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_sigma_dom_realistic))^2 ))/mean(sd_genotypic)*100
sqrt(mean( (sd_genotypic-unlist(sd_ind_dom_realistic))^2 ))/mean(sd_genotypic)*100







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

