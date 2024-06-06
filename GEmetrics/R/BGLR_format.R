#' Format the phenotypic response matrix for BGLR
#'
#' @param Pheno  a data frame with three columns: "Y" for phenotypic values,
#' "Genotype" for genotype names and "Environment" for environment names.
#' All genotypes names must be included in the set or row/column names of
#' the "K" matrix. The number of environments (J) must be at least two
#' @param K a square kinship or genomic relationship matrix for N genotypes
#' whose row and column names include those of the "Genotype" column of
#' the "Pheno" matrix
#' @return a list of two elements: a "BGLR_pheno" phenotypic response matrix
#' with J columns to be used in BGLR and the corresponding "BGLR_K" kinship matrix.
#' @description This function formats the phenotypic data as well as the kinship
#' matrix for BGLR
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
#' ## Generate the phenotypic response matrix for BGLR and the corresponding kinship matrix
#' BGLR_data <- BGLR_format(Pheno=DataSim$Pheno,K=wheat.A)
#' head(BGLR_data$BGLR_pheno)
BGLR_format <- function(Pheno, K){
  ## Number of environments
  J <- length(unique(Pheno$Environment))
  ## Number of genotypes
  N <- nrow(K)
  ## Checks
  if(!all(c("Y","Genotype","Environment")%in%colnames(Pheno))){
    stop("Column names of the Pheno dataframe must include: Y, Genotype and Environment")
  }
  if(J<2L){
    stop("The number of environments must be at least 2")
  }
  if(!any(Pheno$Genotype%in%rownames(K))){
    stop("Some genotype names in the Pheno matrix do not match the row names of the K matrix")
  }
  ## List per environment
  EnvNames <- unique(Pheno$Environment)
  Pheno_list <- lapply(1:J,function(j)Pheno[Pheno$Environment==EnvNames[j],])
  ## Number of replicates per genotype and environment
  rep_mat <- table(factor(Pheno$Genotype,levels = rownames(K)),Pheno$Environment)
  ## BGLR response matrix
  BGLR_pheno <- do.call("rbind",lapply(rownames(K),function(i){
    res <- matrix(NA,nrow = max(rep_mat[i,]),ncol = J,
                  dimnames = list(rep(i,max(rep_mat[i,])),EnvNames))
    G_id <- lapply(1:J,function(j)which(Pheno_list[[j]]$Genotype==i))
    for(j in which(rep_mat[i,]>0)){
      res[1:length(G_id[[j]]),j]=Pheno_list[[j]]$Y[G_id[[j]]]
    }
    return(res)
  }))
  ## BGLR K matrix
  BGLR_K <- K[match(rownames(BGLR_pheno),rownames(K)),match(rownames(BGLR_pheno),colnames(K))]
  return(list(BGLR_pheno=BGLR_pheno,BGLR_K=BGLR_K))
}
