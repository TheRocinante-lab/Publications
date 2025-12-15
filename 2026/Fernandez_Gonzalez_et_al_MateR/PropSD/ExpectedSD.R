SEV2sigma_R <- function(Sigma, muVec, soln) {
  # ---- 0) Minimal checks ----
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma))
    stop("'Sigma' must be a square matrix")
  N <- nrow(Sigma)
  if (length(muVec) != N)
    stop("'muVec' must have length == nrow(Sigma)")
  if (any(soln < 1 | soln > N) || any(soln != floor(soln)))
    stop("'soln' must be integer indices between 1 and N")
  # ---- 1) Unique indices + weights ----
  u    <- unique(soln)                # p unique indices
  wtab <- tabulate(match(soln, u))    # length-p multiplicities
  m    <- length(soln)                # total draws (with repeats)
  # ---- 2) Extract the small p×p subproblem ----
  Sigma_u <- Sigma[u, u, drop = FALSE]  # p×p
  mu_u    <- muVec[u]                   # length-p
  # ---- 3) Compute weighted moments ----
  mu_mean   <- sum(wtab * mu_u) / m
  mu_dev_u  <- mu_u - mu_mean
  
  # ---- 4) SEM^2 via weighted sum of all Sigma_u entries ----
  W_outer  <- wtab %o% wtab          # p×p weight matrix
  SumSigma <- sum(Sigma_u * W_outer)
  SEM2     <- if (SumSigma > 0) SumSigma / (m^2) else 0
  # ---- 5) Row-means in the subproblem (coerce to vector!) ----
  row_means_u <- as.numeric((Sigma_u %*% wtab) / m)  # length-p vector
  # ---- 6) Build Sigma1 on p×p ----
  # Sigma1_ij = Sigma_u[i,j] + SEM2 - row_means_u[i] - row_means_u[j]
  Sigma1 <- Sigma_u +
    SEM2 -
    row_means_u %o% rep(1, length(u)) -
    rep(1, length(u)) %o% row_means_u
  # ---- 7) Compute the two big sums via weighted p×p ops ----
  #  a) sum_s1_sq = sum_{i,j} w_i w_j [Sigma1_ij]^2
  sum_s1_sq <- sum((Sigma1^2) * W_outer)
  #  b) quad form = sum_{i,j} w_i w_j mu_dev_u[i] Sigma1_ij mu_dev_u[j]
  wmudev    <- wtab * mu_dev_u
  quad      <- as.numeric(crossprod(wmudev, Sigma1 %*% wmudev))
  sumAll <- 2 * sum_s1_sq + 4 * quad
  return(sumAll/((m - 1)^2))
}




E_SD <- function(Sigma, soln) {
   
  #The equation for Es2 only works for homogeneous means. We set zero means!
  muVec <- rep(0, nrow(Sigma))
  

  ParentsCount <- sapply(1:nrow(Sigma), function(x) sum(soln == x))
  

  SumDiag <- ParentsCount%*%diag(Sigma) #sum diagonal elements. This is equal to sum(diag(G[allparents,allparents]))
  SumSigma <- t(ParentsCount)%*%Sigma%*%ParentsCount #Sum entire matrix. This is equal to sum(G[allparents,allparents])
  # G_parents <- G[allparents,allparents] #We don't need to compute this huge matrix! Very large memory and time savings!
  SumSigma <- max(1e-6, SumSigma) #this can be lower than zero due to rounding errors, avoid it
  
  #mean of sample variance
  Es2 <- (SumDiag/sum(ParentsCount) - SumSigma/(sum(ParentsCount)^2))*sum(ParentsCount)/(sum(ParentsCount)-1)

  
  #variance of sample variance
  Vars2 <- SEV2sigma_R(Sigma, muVec, soln)

  #Sample variance follows a gamma distribution with the following scale and shape:
  alpha <- Es2^2/Vars2 #shape
  beta <- Vars2/Es2 #scale
  
  return(mean_sqrt_gamma(alpha, beta))
  
}