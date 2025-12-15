library(MateR)

# introduce your MateR license keys here.
# To obtain license keys, please contact j.isidro@upm.es.
# Free license keys will be provided to public bodies,
# such as, Universities, and non-profit organizations.
Username = username
Password = password
Username_TrainSel = username_trainsel
Password_TrainSel = password_trainsel



data(ExampleDataDiploid)



iters <- 1:100


allscenarios <- expand.grid(iters,propsels)
colnames(allscenarios) <- c("iter", "propsel")
dim(allscenarios)

sum(duplicated(paste0(allscenarios[,1],allscenarios[,2]))) #good
unique(allscenarios$iter)#good
table(allscenarios$propsel)#good


#slurm_iter <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) #slurm cluster

#locally:
for (slurm_iter in 1:nrow(allscenarios)) {

  iter <- allscenarios$iter[slurm_iter]
  propsel <- allscenarios$propsel[slurm_iter]
  
  set.seed(1234+iter)
  
  
  
  
  ################################################################################
  #No testers used
  ################################################################################
  #Technow 2012 meffs
  nQTL <- 300
  notQTL <- sample(1:ncol(Markers), (ncol(Markers)-nQTL))
  
  meffs_a <- rgamma(n = ncol(Markers), shape = 0.4, scale = 1.66)
  meffs_a <- meffs_a*sample(c(-1,1), ncol(Markers), replace = T)
  
  meffs_a[notQTL] <- 0
  

  #' 
  #' Dominance degrees were sampled from a normal distribution with mean 
  #' $\mu_{\delta} = 0.2$ and variance $\sigma_{\delta}^2 = 0.2$, representing
  #' typical directional dominance observed in clonally propagated crops 
  #' \citep{wolfe2016genome, endelman2018genetic, yadav2021improved}
  #' 
  #' article{wolfe2016genome,
  #'   title={Genome-wide association and prediction reveals genetic architecture of cassava mosaic disease resistance and prospects for rapid genetic improvement},
  #'   author={Wolfe, Marnin D and Rabbi, Ismail Y and Egesi, Chiedozie and Hamblin, Martha and Kawuki, Robert and Kulakow, Peter and Lozano, Roberto and Carpio, Dunia Pino Del and Ramu, Punna and Jannink, Jean-Luc},
  #'   journal={The Plant Genome},
  #'   volume={9},
  #'   number={2},
  #'   pages={plantgenome2015--11},
  #'   year={2016},
  #'   publisher={Wiley Online Library}
  #' }
  #' @article{endelman2018genetic,
  #'   title={Genetic variance partitioning and genome-wide prediction with allele dosage information in autotetraploid potato},
  #'   author={Endelman, Jeffrey B and Carley, Cari A Schmitz and Bethke, Paul C and Coombs, Joseph J and Clough, Mark E and da Silva, Washington L and De Jong, Walter S and Douches, David S and Frederick, Curtis M and Haynes, Kathleen G and others},
  #'   journal={Genetics},
  #'   volume={209},
  #'   number={1},
  #'   pages={77--87},
  #'   year={2018},
  #'   publisher={Oxford University Press}
  #' }
  #' @article{yadav2021improved,
  #'   title={Improved genomic prediction of clonal performance in sugarcane by exploiting non-additive genetic effects},
  #'   author={Yadav, Seema and Wei, Xianming and Joyce, Priya and Atkin, Felicity and Deomano, Emily and Sun, Yue and Nguyen, Loan T and Ross, Elizabeth M and Cavallaro, Tony and Aitken, Karen S and others},
  #'   journal={Theoretical and Applied Genetics},
  #'   volume={134},
  #'   number={7},
  #'   pages={2235--2252},
  #'   year={2021},
  #'   publisher={Springer}
  #' }
  
  dd_mean <- 0.2 # normal distribution parameter for the dominance degrees
  dd_variance <- 0.2 # normal distribution parameter for the dominance degrees
  
  dd <- rnorm(ncol(Markers), mean = dd_mean, sd = sqrt(dd_variance))
  meffs_d <- abs(meffs_a)*dd
  
  names(meffs_a) <- names(meffs_d) <- colnames(Markers)
  
  # hist(meffs_d[setdiff(1:1000, notQTL)])
  # hist(meffs_a[setdiff(1:1000, notQTL)])
  
  
  
  whichparents <- 1:nrow(Markers)
  size <- 100
  out <- GenomicMatingMT(Parents1 = Parents[whichparents],
                         Parents2 = Parents[whichparents],
                         separator = "/",
                         LD = "Approx",
                         Chromosomes = Chromosomes,
                         ChromosomeCentimorgans = rep(150, length(Chromosomes)),
                         Markers = Markers[whichparents,],
                         parametrization = "Genotypic",
                         phi=4,
                         PropSD = propsel,
                         Markers_T = NULL, #Add testers here
                         markereffects = list(YLD = meffs_a),
                         markereffects_d=list(YLD = meffs_d),
                         n = 0,
                         DHs=F,
                         size = size, #maximum number possible with demo
                         replication = T, #quicker
                         coefficients = NULL,
                         offspring_per_cross = 20,
                         #control = TrainSel::SetControlDefault("demo"),
                         n_selected_per_family=5,
                         Username=Username, #NULL username --> we are using the demo
                         Password=Password, #NULL password --> we are using the demo
                         Username_TrainSel=Username_TrainSel, #NULL username --> we are using the demo
                         Password_TrainSel=Password_TrainSel #NULL password --> we are using the demo
  )
  
  out$OptimalMatingScheme[1:5,]
  out$PropSD
  out$SEM
  out$deltaF
  out$alpha
  out$Optimization_quality
  
  
  Selected <- c()
  for (i in 1:nrow(out$OptimalMatingScheme)) {
    
    Selected <- c(Selected, rep(out$OptimalMatingScheme$Parent1[i],out$OptimalMatingScheme$Number_Of_Crosses[i]),
                  rep(out$OptimalMatingScheme$Parent2[i],out$OptimalMatingScheme$Number_Of_Crosses[i]))
    
  }
  
  parents_AVs <- Markers%*%meffs_a
  
  selected_AVs <- Markers[Selected,]%*%meffs_a
  
  realized_PropSD <- 1-sd(selected_AVs)/sd(parents_AVs)
  
  
  #Genic variance
  Z <- scale(Markers, scale = F)
  Sigma_parents <- t(Z)%*%Z/(nrow(Z)-1)
  sqrt(t(meffs_a) %*% Sigma_parents %*% meffs_a )
  #sd(parents_AVs) #same thing
  
  #genic variance:
  parental_genic_sd <- sqrt(t(meffs_a) %*% diag(diag(Sigma_parents)) %*% meffs_a )
  
  Z2 <- scale(Markers[Selected,], scale = F)
  Sigma_selected <- t(Z2)%*%Z2/(nrow(Z2)-1)
  sqrt(t(meffs_a) %*% Sigma_selected %*% meffs_a )
  #sd(selected_AVs) #same thing
  
  #genic variance:
  selected_genic_sd <- sqrt(t(meffs_a) %*% diag(diag(Sigma_selected)) %*% meffs_a )
  
  realized_PropSD_genic <- 1-selected_genic_sd/parental_genic_sd
  
  
  
  output <- list(Selected = Selected,
                 predicted_PropSD = out$PropSD,
                 realized_PropSD = realized_PropSD,
                 realized_PropSD_genic = realized_PropSD_genic,
                 SEM = out$SEM,
                 deltaF = out$deltaF,
                 alpha = out$alpha)
  output
  #plot(out$TrainSelout$maxvec)
  
  
  save(output, file = paste0("./Results/allResults/Results_line_", iter, "_propsel_", propsel, ".RData"))
  


}


