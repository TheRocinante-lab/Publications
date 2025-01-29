#############################################################################
# Script to optimize the calibration set in genomic selection (maximize the expected reliability).
# Method based on the generalized CD.
# (Isidro et al. 2015)
#############################################################################



###############
#Functions used
###############

# This function creates the matrix of contrast between each of the individual not in the calibration set and the mean of the population
contrasteNonPheno=function(NotSampled)
{
mat=matrix(-1/Nind,Nind,Nind-Nind_in_Sample)
for (i in 1:ncol(mat)) {
mat[NotSampled[i],i]=1-1/Nind
}
return(mat)
}

##############################
# Data required
##########################
matA1=read.table("Amat.csv") #This is the covariance matrix betw the individuals (size Nind x xNind), estimated with the genotypes.
matA1=as.matrix(matA1)
matA2 <- apply(array(colnames(matA1)),1,function(x){strsplit(x,split="X")[[1]][2]})
colnames(matA1) <- matA2
head(matA1)
# validation set of 100
RES=list()
set.seed(123)
il <- 1
#for (il in 1:100){
  
Sample.val<-sample(200,100)
matA1=GBS5[-Sample.val,-Sample.val]

Nind=nrow(matA1) # total number of individuals
nindrep=50 # Choose a size for your calibration set
varP=var(pheno) # Pheno is a vector of phenotypes
h2=0.3	# Trait heritability
varG=h2*varP
varE=(1-h2)/h2*varG
lambda=varE/varG # lambda is needed to estimate the CDmean
lambda=(1-h2)/h2

invA1=solve(matA1) # Inverse of the covariance matrix



##############################
# Optimization algo
##############################

Nind_in_Sample=nindrep
str(Nind_in_Sample)
#Design matrices
Ident<-diag(Nind_in_Sample)
X<-rep(1,Nind_in_Sample)
M<-Ident- (X%*%solve(t(X)%*%X) %*% t(X) )

Sample1<-sample(Nind,Nind_in_Sample) #Calibration set initialization
SaveSample=Sample1
NotSampled1<-seq(1:Nind)
NotSampled<-NotSampled1[-Sample1] # Initial validation set

Z=matrix(0,Nind_in_Sample,Nind)
for (i in 1:length(Sample1)) { Z[i,Sample1[i]]=1 } 

T<-contrasteNonPheno(NotSampled)   # T matrice des contrastes

# Calculate of CDmean of the initial set
matCD<-(t(T)%*%(matA1-lambda*solve(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%matA1%*%T)
CD=diag(matCD)
CDmeanSave=mean(CD)

CDmeanMax1=rep(NA,800)

# Exchange algorithm (maximize CDmean)
cpt2=1
cpt=0
while (cpt2<800) {  # Make sure that 800 is enough in your case (that you reached a plateau), for this look at CDmeanMax1.
 NotSampled=NotSampled1[-Sample1] 
cpt2=cpt2+1
# Remove one individual (randomly choosen) from the sample :
Sample2=sample(Sample1,1)
# Select one individual (randomly choosen) from the individuals that are not in the Calibration set :
Sample3=sample(NotSampled,1)
# New calibration set :
Sample4=c(Sample3,Sample1[Sample1!=Sample2])
# Calculate the mean CD of the new calibration set :
Z=matrix(0,Nind_in_Sample,Nind)
for (i in 1:length(Sample4)) { Z[i,Sample4[i]]=1 } 
NotSampled=NotSampled1[-Sample4] 
T<-contrasteNonPheno(NotSampled)

matCD<-(t(T)%*%(matA1-lambda*solve(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%matA1%*%T)
CD=diag(matCD)

if (mean(CD)>CDmeanSave ) { Sample1=Sample4 # Accept the new Calibration set if CDmean is increased, reject otherwise.
CDmeanSave=mean(CD)  
cpt=0 } else { cpt=cpt+1 
}
CDmeanMax1[cpt2-1]=CDmeanSave
}  #Fin du while

SampleOptimiz=Sample1 # SampleOptimiz is the optimized calibration set

# End

#### for the optimized set repredict the BLUP for each trait
## get the traits in

#tab=read.csv(file="blup07.12.withNA.csv",header=T,as.is=F)

#tab2=tab[match(as.numeric(rownames(gbs.Amat)),tab$GID),-c(1:2)]

### recover correct ID for gbs.Amat considering the removed lines
optim.set=match(rownames(matA1)[SampleOptimiz],rownames(gbs.Amat))

library(rrBLUP)
res=matrix(0,nrow=3,ncol=8)
for (iu in 1:8){
trait=tab2[,iu]
trait2=rep(NA,times=length(trait))
trait2[optim.set]=trait[optim.set]
trait2[Sample.val]=NA
test=mixed.solve(trait2, Z=NULL, K=gbs.Amat, X=NULL, method="REML", 
        bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
        
res[1,iu]=cor(trait[Sample.val],test$u[Sample.val])
### loop 100 time to pick 100 individuals at random and use them as a training population
for (ig in 1:100){
trait2=rep(NA,times=length(trait))
SampleRDm=sample(265,50)

trait2[-Sample.val][SampleRDm]=trait[-Sample.val][SampleRDm]
trait2[Sample.val]=NA
test2=mixed.solve(trait2, Z=NULL, K=gbs.Amat, X=NULL, method="REML", 
        bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
        
res[2,iu]=res[2,iu]+cor(trait[Sample.val],test2$u[Sample.val])
}
res[2,iu]=res[2,iu]/100
### prediction using all the data available
trait=tab2[,iu]
trait2=trait
trait2[Sample.val]=NA
test3=mixed.solve(trait2, Z=NULL, K=gbs.Amat, X=NULL, method="REML", 
        bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
res[3,iu]=cor(trait[Sample.val],test3$u[Sample.val])        
}
RES[[il]]=res
}





