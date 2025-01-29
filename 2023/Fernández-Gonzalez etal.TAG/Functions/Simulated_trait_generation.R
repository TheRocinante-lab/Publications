
simulated_trait_generation <- function (simulated_trait_heritability, geno) {

  simulated_trait_heritability
  markereffects<-rnorm(ncol(geno))
  g<-geno%*%markereffects
  
  #calculate the desired environmental variance from the existing genetic
  #variance and the heritability. This formula is derived from the definition of heritability
  Ve <- as.numeric((var(g)-simulated_trait_heritability*var(g))/simulated_trait_heritability)
  sde <- sqrt(Ve)
  
  #generate the phenotypic data using adding environmental noise to the genetic values
  pheno <- data.frame (GID  = rownames(geno),
                       Trait = g+sde*rnorm(nrow(geno)))
  
  pheno$GID <- factor(as.character(pheno$GID),levels=rownames(K))
  
  heritability_narrow <- as.numeric(var(g)/var(pheno$Trait))
  
  
  cat("Simulated trait heritability:", heritability_narrow, "\n")

  simulated_traits <- list(geno = g, pheno = pheno, h2 = heritability_narrow)
  
  return(simulated_traits)
  
 # save(g, pheno, heritability_narrow, file = filename)
  
    #if the trait is simulated, save the maker effects and the phenotype
    
  
}
