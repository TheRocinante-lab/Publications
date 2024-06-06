#' Compute the best linear unbiased prediction and the conditional variance
#' matrix of environment-specific breeding values
#'
#' @param Pheno a data frame with three columns: "Y" for phenotypic values,
#' "Genotype" for genotype names and "Environment" for environment names.
#' All genotypes names must be included in the set or row and column names of
#' the "K" matrix. The number of environments must be at least two
#' @param K a square kinship or genomic relationship matrix for N genotypes
#' whose row and column names include those of the "Genotype" column of
#' the "Pheno" matrix
#' @param Omega_G a square matrix with genetic covariances between J environments
#' whose row and column names correspond to those of the "Environment" column of
#' the "Pheno" matrix
#' @param Omega_E a square matrix with error covariances between J environments
#' whose row and column names correspond to those of the "Environment" column of
#' the "Pheno" matrix
#' @return a list of two elements: a "G_hat" matrix of best linear unbiased prediction
#' of environment-specific breeding values for the N genotypes (as rows) in J
#' environments (as columns), and a square conditional variance matrix "P"
#' of environment-specific breeding values of dimension NJxNJ where each row
#' block of size N correspond to an environment and the rows of each block
#' correspond to genotypes
#' @description This function calculates the best linear unbiased prediction
#' and the conditional variance matrix of environment-specific breeding values
#' @importFrom stats model.matrix rnorm
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
#' ## Set sparseness by discarding 80% of the combinations
#' Design <- Design[-sample(nrow(Design),round(nrow(Design)*4/5)),]
#'
#' ## Simulate phenotypic data with default parameter values
#' DataSim <- Simulate_MET_data(Design=Design,K=wheat.A)
#'
#' ## Calculate the blup and the conditional variance matrix using simulated variance components
#' ## this step can take several seconds
#' ## note that variance can also be estimated (e.g. using BGLR)
#' BlupEnvBV <- EnvBV_blup(Pheno=DataSim$Pheno,K=wheat.A,Omega_G=DataSim$Omega_G,
#'                         Omega_E=DataSim$Omega_E)
#'
#' ## Display results
#' head(BlupEnvBV$G_hat)
#' BlupEnvBV$P[1:5,1:5]
#'
EnvBV_blup <- function(Pheno, K, Omega_G, Omega_E){
  ## Number of environments
  J <- ncol(Omega_G)
  ## Number of genotypes
  N <- nrow(K)
  ## Checks
  if(!all(c("Y","Genotype","Environment")%in%colnames(Pheno))){
    stop("Column names of the Pheno dataframe must include: Y, Genotype and Environment")
  }
  if(length(unique(Pheno$Environment))<2L){
    stop("The number of environments must be at least 2")
  }
  if(!any(Pheno$Genotype%in%rownames(K))){
    stop("Some genotype names in the Pheno matrix do not match the row names of the K matrix")
  }
  if(nrow(Omega_G)!=J|ncol(Omega_G)!=J|
     !all(sort(rownames(Omega_G))==sort(unique(Pheno$Environment)))|
     !all(sort(colnames(Omega_G))==sort(unique(Pheno$Environment)))){
    stop("Some environment names do not match between the Environment column of the Pheno data frame and the row/column names of the Omega_G matrix")
  }
  if(nrow(Omega_E)!=J|ncol(Omega_E)!=J|
     !all(sort(rownames(Omega_E))==sort(unique(Pheno$Environment)))|
     !all(sort(colnames(Omega_E))==sort(unique(Pheno$Environment)))){
    stop("Some environment names do not match between the  Environment column of the Pheno data frame and the row/column names of the Omega_E matrix")
  }
  ## Covariance matrices
  Sigma_G <- kronecker(Omega_G,K,make.dimnames = T)
  Sigma_E <- diag(nrow(Pheno))
  EnvNames <- unique(Pheno$Environment)
  for(j in 1:J){
    for(jj in 1:J){
      Sigma_E_j_jj <- Sigma_E[Pheno$Environment==EnvNames[j],Pheno$Environment==EnvNames[jj]]
      rownames(Sigma_E_j_jj) <- Pheno$Genotype[Pheno$Environment==EnvNames[j]]
      colnames(Sigma_E_j_jj) <- Pheno$Genotype[Pheno$Environment==EnvNames[jj]]
      for(i in unique(Pheno$Genotype)){
        if(i%in%rownames(Sigma_E_j_jj) & i%in%colnames(Sigma_E_j_jj)){
          Sigma_E_j_jj[i,i] <- Omega_E[j,jj]
        }
      }
    }
  }
  ## Design matrices
  X <- stats::model.matrix(~factor(Pheno$Environment,levels = EnvNames)-1)
  Z <- stats::model.matrix(~factor(apply(Pheno[,c("Environment","Genotype")],1,paste,collapse=":"),
                            levels = rownames(Sigma_G))-1)
  ## Phenotypic covariance matrix
  Sigma_Y <- Z%*%Sigma_G%*%t(Z) + Sigma_E
  Sigma_Y_inv <- solve(Sigma_Y)
  ## Best linear unbiased estimates
  beta_hat <- solve(t(X)%*%Sigma_Y_inv%*%X)%*%t(X)%*%Sigma_Y_inv%*%Pheno$Y
  ## Conditional covariance matrix
  P <- Sigma_G - Sigma_G%*%t(Z)%*%Sigma_Y_inv%*%Z%*%Sigma_G
  ## Best linear unbiased prediction
  M <- Sigma_Y_inv - Sigma_Y_inv%*%X%*%solve(t(X)%*%Sigma_Y_inv%*%X)%*%t(X)%*%Sigma_Y_inv
  G_hat <- matrix(Sigma_G%*%t(Z)%*%M%*%Pheno$Y,N,J,dimnames = list(rownames(K),EnvNames)) +
    tcrossprod(rep(1,N),beta_hat)
  return(list(G_hat=G_hat,P=P))
}

