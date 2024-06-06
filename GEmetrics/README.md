
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEmetrics

<!-- badges: start -->
<!-- badges: end -->

The goal of GEmetrics is to provide functions to calculate the best
linear unbiased prediction (BLUP) of the following
genotype-by-environment (GE) metrics: ecovalence, environmental
variance, Finlay and Wilkinson regression and Lin and Binns superiority
measure, based on a multi-environment genomic prediction model.

## Installation

You can install GEmetrics directly from the CRAN:

``` r
install.packages(pkg='GEmetrics',repos='https://cran.r-project.org/')
```

or from GitHub:

``` r
install.packages(pkg='devtools',repos='https://cran.r-project.org/')     ## install devtools
devtools::install_git('https://github.com/TheRocinante-lab/GEmetrics')   ## install GEmetrics from GitHub
```

## Example

This is a basic example which shows you how to calculate the BLUP of GE
metrics.

#### Simulate phenotypic data

Multi-environment trial data is first simulated based on the “wheat”
dataset from BGLR.

A design data frame is generated displaying all combinations of
genotypes and environments. Some combinations are discarded to set
sparseness in the data (75% here).

``` r
## Set seed for reproductibility
set.seed(123)

## Load "wheat" dataset from BGLR
data("wheat",package = "BGLR")

## Generate a design data frame for all genotypes in 5 environments
Design <- expand.grid(Genotype=rownames(wheat.A),Environment=paste0("Env",1:5))

## Set sparseness by discarding 75% of the combinatons
Design <- Design[-sample(nrow(Design),round(nrow(Design)*3/4)),]
head(Design)
#>    Genotype Environment
#> 3      2167        Env1
#> 6      3889        Env1
#> 15    13396        Env1
#> 22    14103        Env1
#> 26    16004        Env1
#> 28    16262        Env1
```

Phenotypes are then simulated using trait and environments parameters:

- h2: heritability, either a single value for a heritability common to
  all environments (e.g. 0.5 in the example below), or a vector of
  heritabilities associated with each environment
- rho: genetic correlations between environment pairs, either a single
  value for a genetic correlation common to all environment pairs
  (e.g. 0.5 in the example below), or a square correlation matrix
- sd_mu: standard deviation of the Gaussian distribution in which
  environment means are drawn (e.g. 1 in the example below)

``` r
## Simulate phenotypic data with default parameter values
DataSim <- GEmetrics::Simulate_MET_data(Design=Design,K=wheat.A,h2=0.5,rho=0.5,sd_mu=1)
```

The resulting DataSim object include:

1)  Pheno: data frame with simulated phenotypes

``` r
## Simulated phenotypes
head(DataSim$Pheno)
#>             Y Genotype Environment
#> 3   1.5105722     2167        Env1
#> 6   0.1510397     3889        Env1
#> 15 -3.3939121    13396        Env1
#> 22  0.3920105    14103        Env1
#> 26 -0.4409719    16004        Env1
#> 28 -1.1378441    16262        Env1
```

2)  EnvBV: matrix of simulated environment-specific breeding values

``` r
## Simulated environment-specific breeding values
head(DataSim$EnvBV)
#>            Env1       Env2       Env3       Env4      Env5
#> 775  -0.8188854  0.7310795 -0.7722164 -2.9188405 0.1487880
#> 2166  0.5865631  2.1463691  1.9910743  0.7155801 1.7982954
#> 2167  0.6615584  2.1808788  1.9865398  0.7955096 1.9156231
#> 2465  0.6525278  0.7716421  0.5629410  0.0811664 1.8413287
#> 3881  0.3101735  0.6611265  1.6816084 -2.2797745 1.7717076
#> 3889 -1.0755458 -2.0130535 -1.5269430 -2.4387551 0.5240766
```

3)  Omega_G: genetic covariance matrix between environments

``` r
## Genetic covariance matrix between environments
DataSim$Omega_G
#>      Env1 Env2 Env3 Env4 Env5
#> Env1  1.0  0.5  0.5  0.5  0.5
#> Env2  0.5  1.0  0.5  0.5  0.5
#> Env3  0.5  0.5  1.0  0.5  0.5
#> Env4  0.5  0.5  0.5  1.0  0.5
#> Env5  0.5  0.5  0.5  0.5  1.0
```

4)  Omega_E: error covariance matrix between environments

``` r
## Error covariance matrix between environments
DataSim$Omega_E
#>      Env1 Env2 Env3 Env4 Env5
#> Env1    1    0    0    0    0
#> Env2    0    1    0    0    0
#> Env3    0    0    1    0    0
#> Env4    0    0    0    1    0
#> Env5    0    0    0    0    1
```

#### Estimate variance components using BGLR

From simulated data, variance components can be estimated using an
inference method like [BGLR](https://github.com/gdlc/BGLR-R), or any
other methods able to infer Omega_G and Omega_E.

First, the phenotypic data frame must be transformed into a phenotypic
response matrix.

``` r
## Generate the phenotypic response matrix for BGLR and the corresponding K matrix
BGLR_data <- GEmetrics::BGLR_format(Pheno=DataSim$Pheno,K=wheat.A)
head(BGLR_data$BGLR_pheno)
#>           Env1        Env2        Env3      Env4     Env5
#> 775         NA -0.04697659          NA        NA       NA
#> 2167 1.5105722          NA          NA  1.580189 1.205947
#> 2465        NA          NA  1.49828160        NA       NA
#> 3881        NA  0.50301295          NA        NA       NA
#> 3889 0.1510397          NA          NA        NA       NA
#> 4248        NA          NA -0.04786517 -1.401185       NA
```

The inference can be done using the “Multitrait” function of BGLR to
estimate the Omega_G and Omega_E covariance matrices. Note that the
current CRAN version of BGLR (October 2023) may lead to an issue when
the phenotypic data is very sparse, but not the most recent GitHub
version.

``` r
## Run BGLR inference
ETA<-list(list(K=BGLR_data$BGLR_K,model="RKHS"))
BGLR_results <- BGLR::Multitrait(y=BGLR_data$BGLR_pheno,ETA=ETA,
                                 resCov=list(type="DIAG"),
                                 nIter=1000,burnIn=500,verbose = F,saveAt = "Test_")
#> Checking variance co-variance matrix K  for linear term 1
#> Ok
#> Setting linear term 1
#> MSx=1.60381006430079
#> UNstructured covariance matrix
#> df0 was set to 6
#> S0 set to
#>            Env1      Env2     Env3     Env4      Env5
#> Env1 10.1237857  3.523954 2.648470 2.804323 0.9866993
#> Env2  3.5239541 11.076081 2.746944 3.647022 1.7981729
#> Env3  2.6484702  2.746944 9.101782 2.333046 2.4699484
#> Env4  2.8043233  3.647022 2.333046 9.662190 1.1637173
#> Env5  0.9866993  1.798173 2.469948 1.163717 5.9543964
#> Initializing resCov
#> Setting hyperparameters for DIAG R
#> df0 set to  5 for all the traits
#> S0 was set to
#>      Env1      Env2      Env3      Env4      Env5 
#>  9.471367 10.362293  8.515226  9.039518  5.570671
#> Done
unlink(c("Test_R.dat","Test_Omega_1.dat","Test_mu.dat"))
Omega_G <- BGLR_results$ETA[[1]]$Cov$Omega
Omega_E <- BGLR_results$resCov$R
rownames(Omega_E) <- rownames(Omega_E) <- rownames(Omega_G)
```

The estimate of the genetic covariance matrix Omega_G is:

``` r
Omega_G
#>           Env1      Env2      Env3      Env4      Env5
#> Env1 1.0247664 0.4677599 0.3996223 0.4624807 0.2011146
#> Env2 0.4677599 0.9498619 0.3683848 0.3746494 0.2067460
#> Env3 0.3996223 0.3683848 0.7740551 0.3278520 0.2372795
#> Env4 0.4624807 0.3746494 0.3278520 0.8901805 0.2105031
#> Env5 0.2011146 0.2067460 0.2372795 0.2105031 0.4848253
```

and the estimate of the error covariance matrix Omega_E is:

``` r
Omega_E
#>          [,1]     [,2]     [,3]    [,4]      [,5]
#> Env1 1.118006 0.000000 0.000000 0.00000 0.0000000
#> Env2 0.000000 1.045441 0.000000 0.00000 0.0000000
#> Env3 0.000000 0.000000 1.141221 0.00000 0.0000000
#> Env4 0.000000 0.000000 0.000000 1.21129 0.0000000
#> Env5 0.000000 0.000000 0.000000 0.00000 0.9782482
```

#### Calculate BLUP and conditional variance of environment-specific breeding values

The BLUP and the conditional variance of environment-specific breeding
values can be calculated from the phenotypes and the variance component
estimates. Note that the BLUPs could also be obtained directly from BGLR
outputs.

``` r
## Calculate BLUP and conditional variance
BlupEnvBV <- GEmetrics::EnvBV_blup(Pheno=DataSim$Pheno,K=wheat.A,Omega_G=Omega_G,Omega_E=Omega_E)
```

The BLUPs obtained:

``` r
head(BlupEnvBV$G_hat)
#>             Env1       Env2       Env3       Env4      Env5
#> 775  -0.66753201 -0.2112071  0.4905182 -1.5057121 0.8578696
#> 2166  0.81771762  0.4090425  1.3872625  0.5760052 1.4273601
#> 2167  0.82168849  0.4112842  1.3889899  0.5801013 1.4281564
#> 2465  0.64992740  0.9668607  1.4473134 -0.3386963 1.5218385
#> 3881 -0.08437372  0.2692015  0.9020293 -1.6997681 0.6887119
#> 3889 -1.37042995 -2.3024314 -0.5541623 -2.0430040 0.4299118
```

and the conditional variance matrix:

``` r
BlupEnvBV$P[1:5,1:5]
#>             Env1:775    Env1:2166    Env1:2167    Env1:2465    Env1:3881
#> Env1:775  1.09832586  0.013557905  0.013513339  0.053213436  0.241726662
#> Env1:2166 0.01355791  0.517089943  0.514962943 -0.001668494 -0.006067685
#> Env1:2167 0.01351334  0.514962943  0.516518608 -0.001649603 -0.006042546
#> Env1:2465 0.05321344 -0.001668494 -0.001649603  0.502314319  0.015195531
#> Env1:3881 0.24172666 -0.006067685 -0.006042546  0.015195531  0.776973074
```

#### Obtain GE metrics estimates

Each GE metric can be estimated using the complete BLUP involving both
the squared expectation and the variance term:

``` r
metrics <- c("Ecovalence","EnvironmentalVar","FinlayWilkRegression","LinBinns")
GEmetrics_hat_geno_exp_var <- sapply(metrics,function(m)
  GEmetrics::GEmetrics_blup(G_hat=BlupEnvBV$G_hat,metric=m,P=BlupEnvBV$P))
head(GEmetrics_hat_geno_exp_var)
#>      Ecovalence EnvironmentalVar FinlayWilkRegression LinBinns
#> 775    2.394463        1.4653910            1.0064156 5.107793
#> 2166   3.337781        0.7363370            0.4428062 2.576045
#> 2167   3.347547        0.7364124            0.4414242 2.572119
#> 2465   1.631830        0.9115902            0.7942961 2.281372
#> 3881   2.532864        1.5401010            1.0298445 4.236280
#> 3889   1.780832        1.4125020            1.0651308 8.526201
```

or using the partial BLUP including the the squared expectation only:

``` r
metrics <- c("Ecovalence","EnvironmentalVar","FinlayWilkRegression","LinBinns")
GEmetrics_hat_geno_exp <- sapply(metrics,function(m)
  GEmetrics::GEmetrics_blup(G_hat=BlupEnvBV$G_hat,metric=m))
head(GEmetrics_hat_geno_exp)
#>      Ecovalence EnvironmentalVar FinlayWilkRegression LinBinns
#> 775   0.1037494        0.8798098            1.0043964 4.383918
#> 2166  1.3178871        0.2163899            0.4332007 1.927741
#> 2167  1.3237323        0.2154950            0.4318088 1.923486
#> 2465  0.3207806        0.5687083            0.7885727 1.797779
#> 3881  0.6793336        1.0639080            1.0281445 3.616348
#> 3889  1.2136119        1.2559158            1.0626649 8.104803
```

Estimated can be compared to the true GE metric values obtained from
simulated environment-specific breeding values using the correlation:

``` r
GEmetrics_true <- sapply(metrics,function(m)GEmetrics::GEmetrics_blup(G_hat=DataSim$EnvBV,metric=m,P=NULL))
data.frame("Geno_Exp_Var" = diag(cor(GEmetrics_hat_geno_exp_var,GEmetrics_true)),
           "Geno_Exp" = diag(cor(GEmetrics_hat_geno_exp,GEmetrics_true)))
#>                      Geno_Exp_Var  Geno_Exp
#> Ecovalence              0.3596291 0.2137034
#> EnvironmentalVar        0.5047009 0.4749053
#> FinlayWilkRegression    0.5315392 0.5316567
#> LinBinns                0.7944993 0.7935010
```

or the root mean-square error of estimation:

``` r
data.frame("Geno_Exp_Var" = sqrt(colMeans((GEmetrics_hat_geno_exp_var-GEmetrics_true)^2)),
           "Geno_Exp" = sqrt(colMeans((GEmetrics_hat_geno_exp-GEmetrics_true)^2)))
#>                      Geno_Exp_Var  Geno_Exp
#> Ecovalence              2.0265446 3.0840888
#> EnvironmentalVar        0.8941891 1.1087231
#> FinlayWilkRegression    0.3837222 0.3837641
#> LinBinns                3.4556773 4.0040074
```
