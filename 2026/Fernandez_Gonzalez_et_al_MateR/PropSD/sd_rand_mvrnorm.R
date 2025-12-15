library(MateR)
library(MASS)
source("ExpectedSD.R")

#line
data(ExampleDataDiploid)


Sigma <- G

n_reps <- 1e5

iters <- 1:20
solnSizes <- c(5,20,100,200)
scenarios <- expand.grid(iters,solnSizes)
colnames(scenarios) <- c("iters", "sizes")
dim(scenarios)

outdf <- data.frame()

set.seed(1)
for (i in 1:nrow(scenarios)) {
  
  print(i)
  size <- scenarios$sizes[i]
  soln <- sample(1:nrow(Sigma), size, replace = T)
  
  

  SampleWithCorrelations2 <- mvrnorm(n_reps, mu = rep(0,nrow(Sigma)), Sigma = Sigma)
  SampleWithCorrelations <- SampleWithCorrelations2[,soln]
  sdSelected <- apply(SampleWithCorrelations, 1, function(x) {
    return(sd(x))
  }) 
  sdAll <- apply(SampleWithCorrelations2, 1, function(x) {
    return(sd(x))
  }) 
  
  propRealized <- 1-sdSelected/sdAll
  
  
  
  propExpected <- 1-E_SD(Sigma, soln)/E_SD(Sigma, 1:nrow(Sigma))
  
  
  
  
  tmp <- data.frame(scenario = i,
                    size = size,
                    propRealized = propRealized,
                    propExpected = propExpected)
  
  
  

  outdf <- rbind(outdf, tmp)
  
}

colnames(outdf) <- c("scenario", "size", "propRealized", "propExpected")



save(outdf, file = "./Results/Diploid_ResultsRand.RData")








#clonal
rm(list=ls())
data(ExampleDataAutotetraploid)


Sigma <- G4

n_reps <- 1e5

iters <- 1:20
solnSizes <- c(5,20,100,200)
scenarios <- expand.grid(iters,solnSizes)
colnames(scenarios) <- c("iters", "sizes")


outdf4 <- data.frame()

set.seed(1234)
for (i in 1:nrow(scenarios)) {
  
  print(i)
  size <- scenarios$sizes[i]
  soln <- sample(1:nrow(Sigma), size, replace = T)
  
  
  
  SampleWithCorrelations2 <- mvrnorm(n_reps, mu = rep(0,nrow(Sigma)), Sigma = Sigma)
  SampleWithCorrelations <- SampleWithCorrelations2[,soln]
  sdSelected <- apply(SampleWithCorrelations, 1, function(x) {
    return(sd(x))
  }) 
  sdAll <- apply(SampleWithCorrelations2, 1, function(x) {
    return(sd(x))
  }) 
  
  propRealized <- 1-sdSelected/sdAll
  
  
  
  propExpected <- 1-E_SD(Sigma, soln)/E_SD(Sigma, 1:nrow(Sigma))
  
  
  
  
  tmp <- data.frame(scenario = i,
                    size = size,
                    propRealized = propRealized,
                    propExpected = propExpected)
  
  

  
  outdf4 <- rbind(outdf4, tmp)
  
}

colnames(outdf4) <- c("scenario", "size", "propRealized", "propExpected")


save(outdf4, file = "./Autotetraploid_ResultsRand.RData")


