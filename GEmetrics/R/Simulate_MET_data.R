#' Simulate multi-environment trials data
#'
#' @param Design a data frame with two columns: "Genotype" for genotype names and
#' "Environment" for environment names. All genotypes names must be included
#' in the set of row/column names of the "K" matrix. The number of environments (J)
#' must be at least two
#' @param K a square kinship or genomic relationship matrix for N genotypes
#' whose row/column names include those of the "Genotype" column of
#' the "Design" matrix
#' @param h2 heritability (numeric value between 0 and 1 excluded) of observations
#' in each environment: either a scalar to set a common heritability for all
#' environments, or a vector of heritabilities associated with each environment
#' of size J. The default value is a heritability of 0.5 for all environments
#' @param rho genetic correlation (numeric value between -1 and 1 excluded)
#' between environment pairs: either a scalar to set a common genetic correlation
#' between all environment pairs, or a square correlation matrix of dimension JxJ.
#' The default value is a genetic correlation of 0.5 between all environment pairs
#' @param sd_mu standard deviation (positive numeric value) of the Gaussian
#' distribution in which environment means are drawn. The default value is 1
#' @return a list of two elements: a "Pheno" data frame consisting of the
#' "Design" data frame to which a "Y" column containing simulated phenotypic
#'  values has been added, and a "EnvBV" matrix of dimension NxJ containing the
#'  simulated environment-specific breeding values
#' @description This function calculates the best linear unbiased prediction
#' and the conditional variance matrix of environment-specific breeding values
#' @importFrom BGLR BGLR
#' @export
#'
#' @examples
#' ## Set seed for reproductibility
#' set.seed(123)
#'
#' ## Load "wheat" dataset from BGLR
#' data("wheat",package = "BGLR")
#'
#' ## Generate a design data frame for all genotypes in 5 environments
#' Design <- expand.grid(Genotype=rownames(wheat.A),Environment=paste0("Env",1:5))
#'
#' ## Set sparseness by discarding 75% of the combinations
#' Design <- Design[-sample(nrow(Design),round(nrow(Design)*3/4)),]
#'
#' ## Simulate phenotypic data with default parameter values
#' DataSim <- Simulate_MET_data(Design=Design,K=wheat.A)
#'
#' ## Simulated phenotypes
#' head(DataSim$Pheno)
#'
#' ## Simulated environment-specific breeding values
#' head(DataSim$EnvBV)
#'
#' ## Genetic covariance matrix between environments
#' DataSim$Omega_G
#'
#' ## Error covariance matrix between environments
#' DataSim$Omega_E
#'
Simulate_MET_data <- function(Design,K,h2=0.5,rho=0.5,sd_mu=1){
  ## Number of genotypes
  N <-  nrow(K)
  ## Number of environments
  J <- length(unique(Design$Environment))
  ## Checks
  if(!all(c("Genotype","Environment")%in%colnames(Design))){
    stop("Column names of the Design dataframe must include: Genotype and Environment")
  }
  if(J<2L){
    stop("The number of environments must be at least 2")
  }
  if(!is.double(rho)|(is.vector(rho)&length(rho)>1)){
    stop("The rho parameter must be a numeric value specified as a scalar or square matrix")
  }else if(is.matrix(rho)){
    if(length(rho)!=J^2){
      stop("When specifying a matrix of values for rho, the matrix dimension must match the number of environments in the Design matrix")
    }else if(!isSymmetric.matrix(rho)){
      stop("When specifying a matrix of values for rho, the matrix must be symmetric")
    }else if(any(diag(rho)!=1)){
      stop("When specifying a matrix of values for rho, the diagonal must only include 1")
    }else if(any(eigen(rho)$values<=0)){
      stop("When specifying a matrix of values for rho, it must be positive definite")
    }
  }else{
    if(rho>=1|rho<=(-1)){
      stop("The rho parameter must be between -1 and 1 excluded")
    }else if(any(eigen(matrix(rho,J,J)+diag(J)*(1-rho))$values<=0)){
      stop("When specifying a common rho between all environment pairs, the resulting genetic correlation matrix must be positive definite")
    }
  }
  if(!is.double(h2)|is.matrix(h2)){
    stop("The h2 parameter must be a numeric value specified as a scalar or vector")
  }else if(any(h2<=0)|any(h2>=1)){
    stop("The h2 parameter must be between 0 and 1 excluded")
  }else if(length(h2)>1L&length(h2)!=J){
      stop("When specifying a vector of values for h2, the vector size must match the number of environments in the Design matrix")
  }
  if(!is.double(sd_mu)|length(sd_mu)!=1){
    stop("The sd_mu parameter must be a numeric value specified as a scalar")
  }else if(sd_mu<0){
    stop("The sd_mu parameter must be positive")
  }
  if(!any(Design$Genotype%in%rownames(K))){
    stop("Some genotype names in the Design matrix do not match the row names of the K matrix")
  }
  ## Cholesky factorization of the K matrix
  chol_K <- chol(K+diag(N)*0.001)
  ## Environment means
  mu_env <- rnorm(J,0,sd_mu)
  ## Genetic correlations
  if(!is.matrix(rho)){
    rho <- matrix(rho,J,J) + diag(J)*(1-rho)
  }
  ##Heritability
  if(length(h2)==1){
    h2 <- rep(h2,J)
  }
  ## Genetic covariance matrix between environments
  Omega_G <- rho
  rownames(Omega_G) <- unique(Design$Environment)
  colnames(Omega_G) <- unique(Design$Environment)
  ## Error covariance matrix between environments
  Omega_E <- diag((1-h2)/h2)
  rownames(Omega_E) <- unique(Design$Environment)
  colnames(Omega_E) <- unique(Design$Environment)
  ## Env-BVs
  EnvBV <- matrix(rnorm(N*J),J,N)%*%chol_K
  EnvBV <- tcrossprod(rep(1,N),mu_env) + t(EnvBV)%*%chol(rho)%*%diag(rep(1,J))
  rownames(EnvBV) <- rownames(K)
  colnames(EnvBV) <- unique(Design$Environment)
  ## Design matrix
  Z <- model.matrix(~factor(Design$Genotype,levels=rownames(K))-1)
  G_init <- Z%*%EnvBV
  G <- sapply(1:nrow(Design),function(i)
    G_init[i,which(unique(Design$Environment)%in%Design$Environment[i])])
  ## Errors
  E_list <- lapply(unique(Design$Environment),function(j){
    j_size <- length(which(Design$Environment==j))
    E_j <- rnorm(j_size,0,sqrt(diag(Omega_E)))
    return(E_j)
  })
  E <- rep(0,nrow(Design))
  for(j in 1:J){
    E[Design$Environment==unique(Design$Environment)[j]] = E_list[[j]]
  }
  ## Phenotypes
  Y <- G + E
  Pheno <- data.frame(Y = Y, Design[,c("Genotype","Environment")])
  return(list(Pheno=Pheno, EnvBV=EnvBV,Omega_G=Omega_G,Omega_E=Omega_E))
}

