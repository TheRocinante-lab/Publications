#' Best linear unbiased prediction of genotype-by-environment (GE) metrics
#'
#' @param G_hat a matrix of best linear unbiased prediction of environment-specific
#' breeding values for N genotypes (as rows) in J environments (as columns)
#' @param metric a character string indicating what GE metric to consider:
#' "Ecovalence", "EnvironmentalVar", "FinlayWilkRegression", "LinBinns"
#' @param P (optional) a square conditional variance matrix of environment-specific
#' breeding values of dimension NJxNJ, where each row block of size N correspond
#' to an environment and the rows of each block correspond to genotypes
#' @return a vector of size N with best linear unbiased prediction of
#' the GE metric
#' @description This function calculates the best linear unbiased prediction of
#' the following GE metrics: ecovalence, environmental variance, Finlay and
#' Wilkinson regression and Lin and Binns superiority measure.
#' Ignoring the P matrix resumes to ignoring the condition
#' variance term in the calculation
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
#'
#' ## Compute GE metric estimates
#' GEmetrics_hat <- GEmetrics_blup(G_hat=BlupEnvBV$G_hat,metric="Ecovalence",P=BlupEnvBV$P)
#' head(GEmetrics_hat)
GEmetrics_blup <- function(G_hat, metric, P=NULL){
  ## Number of environments
  J <- ncol(G_hat)
  ## Number of genotypes
  N <- nrow(G_hat)
  ## Checks
  if(!metric%in%c("Ecovalence", "EnvironmentalVar", "FinlayWilkRegression", "LinBinns")){
    stop("The metric field must be one of the following: Ecovalence, EnvironmentalVar, FinlayWilkRegression, LinBinns")
  }
  if(!is.null(P)){
    if(nrow(P)!=J*N|ncol(P)!=J*N){
      stop("The P matrix must be a square matrix of dimension NJxNJ where N and J are the number of genotypes and environments in the G_hat matrix, respectively")
    }
  }
  if(!is.numeric(G_hat)|!is.matrix(G_hat)){
    stop("G_hat must be a numeric matrix of dimension NxJ where N and J are the number of genotypes and environments in the G_hat matrix, respectively")
  }
  ## Test which metric to calculate
  if(metric=="Ecovalence"){
    ## Squared expectation term
    res <- rowSums((G_hat-rowMeans(G_hat)-tcrossprod(rep(1,N),colMeans(G_hat))+mean(G_hat))^2)
    ## Conditional variance term
    if(!is.null(P)){
      V_Y1a <- rowSums(matrix(diag(P),ncol = J))
      P_bdiag <- lapply(1:J,function(t)P[(1+N*(t-1)):(N*t),(1+N*(t-1)):(N*t)])
      V_Y1b <- rowSums(sapply(P_bdiag,rowSums))
      V_Y1c <- sum(sapply(P_bdiag,sum))
      P_b <- lapply(1:J,function(i)lapply(1:J,function(e)P[(1+N*(i-1)):(N*i),(1+N*(e-1)):(N*e)]))
      V_Y2a <- rowSums(sapply(P_b,function(i)rowSums(sapply(i,diag))))
      V_Y2b <- rowSums(sapply(P_b,function(i)rowSums(sapply(i,rowSums))))
      V_Y2c <- sum(sapply(P_bdiag,function(i)sum(sapply(i,sum))))
      res <- res + (V_Y1a -2/N*V_Y1b + 1/(N^2)*V_Y1c) - 1/J*(V_Y2a -2/N*V_Y2b + 1/(N^2)*V_Y2c)
    }
  }else if(metric=="EnvironmentalVar"){
    ## Squared expectation term
    res <- 1/(J-1)*rowSums((G_hat-rowMeans(G_hat))^2)
    ## Conditional variance term
    if(!is.null(P)){
      res_V1 <- 1/(J-1)*rowSums(matrix(diag(P),ncol = J))
      res_V2 <- 1/((J-1)*J)*rowSums(Reduce("+",lapply(1:J,function(i)sapply(1:J,function(e)
        diag(P[(1+N*(i-1)):(N*i),(1+N*(e-1)):(N*e)])))))
      res <- res + res_V1 - res_V2
    }
  }else if(metric=="FinlayWilkRegression"){
    ## Squared expectation term
    if(is.null(P)){
      res <- rowSums((G_hat-rowMeans(G_hat))*tcrossprod(rep(1,N),(colMeans(G_hat)-mean(G_hat))))/
        sum((colMeans(G_hat)-mean(G_hat))^2)
    }else{
      ## Conditional variance term
      Denom_E2_Y <- sum((colMeans(G_hat)-mean(G_hat))^2)
      P_bdiag <- lapply(1:J,function(t)P[(1+N*(t-1)):(N*t),(1+N*(t-1)):(N*t)])
      P_b <- lapply(1:J,function(i)lapply(1:J,function(e)P[(1+N*(i-1)):(N*i),(1+N*(e-1)):(N*e)]))
      Denom_V_Y1 <- sum(unlist(P_bdiag))
      Denom_V_Y2 <- sum(P)
      Denom <- Denom_E2_Y + 1/(N^2)*Denom_V_Y1 - 1/(N^2*J)*Denom_V_Y2
      Num_E2_Y <- rowSums((G_hat-rowMeans(G_hat))*tcrossprod(rep(1,N),(colMeans(G_hat)-mean(G_hat))))
      Num_Cov_Y1 <- rowSums(sapply(P_bdiag,rowSums))
      Num_Cov_Y2 <- rowSums(sapply(P_b,function(i)rowSums(sapply(i,rowSums))))
      Num <- Num_E2_Y + 1/N*Num_Cov_Y1 - 1/(N*J)*Num_Cov_Y2
      res <- Num/Denom
    }
  }else if(metric=="LinBinns"){
    ## Squared expectation term
    ref <- apply(G_hat,2,which.max)
    res <- 1/(J*2)*rowSums((G_hat-tcrossprod(rep(1,N),diag(G_hat[ref,])))^2)
    ## Conditional variance term
    if(!is.null(P)){
      P_list <- lapply(1:N,function(i)lapply(1:J,function(j)
        P[c(i+N*(j-1),ref[j]+N*(j-1)),c(i+N*(j-1),ref[j]+N*(j-1))]))
      res <- res + 1/(J*2)*sapply(P_list,function(i)sum(sapply(i,function(j)
        sum(2*diag(diag(j))-j))))
    }
  }

  return(res)
}

