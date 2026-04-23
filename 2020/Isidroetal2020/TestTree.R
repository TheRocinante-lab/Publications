depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  } else{
    return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))    
  }
}



remidpoint <- function(x){
  isdendro <- inherits(x, "dendrogram")
  setnodeattr <- function(node){
    if(is.list(node)){
      get_cladesizes <- function(z) length(unlist(z, use.names = FALSE))
      cladesizes <- vapply(node, get_cladesizes, 0)
      nclades <- length(cladesizes)
      attr(node, "members") <- sum(cladesizes)
      attr(node, "midpoint") <- ((cladesizes[1] - 1)/2 +
                                   #(cladesizes[1] +
                                   (sum(cladesizes[1:(nclades - 1)]) +
                                      (cladesizes[nclades] - 1)/2))/2
      attr(node, "leaf") <- NULL
    }else{
      attr(node, "members") <- 1
      attr(node, "leaf") <- TRUE
    }
    return(node)
  }
  settreeattr <- function(tree){
    tree <- setnodeattr(tree)
    if(is.list(tree)) tree[] <- lapply(tree, settreeattr)
    return(tree)
  }
  x <- settreeattr(x)
  if(isdendro) class(x) <- "dendrogram"
  return(x)
}



TestTree<-function(dend, SNPs,Pheno, nlevels=5,alpha=.05, conditional=TRUE, correction="PC", nPC=10){
  
  
  nodefunc<-function(n){if (is.leaf(n)){
    return(labels(n))
  } else {NA}
  }
  outerfuncmemberlabels <- function(n) {
    a <- attributes(n)
    attr(n, "nodePar") <-
      c(a$nodePar, list(memberlabels=c(unlist(dendrapply(n, nodefunc)))))
    n
  }
  
  depthtree<-depth(dend)
  if (nlevels>depthtree){nlevels<-depthtree}
  dend <- dendrapply(dend, outerfuncmemberlabels)
  
  allmarkers<-colnames(SNPs)
  
  outerfuncoutermarkers <- function(n) {
    outmarkers<-setdiff(allmarkers,attr(n, "nodePar")$memberlabels)
    attr(n, "nodePar") <-
      c(attr(n, "nodePar"), list(outmarkers=outmarkers))
    n
  }
  
  dend <- dendrapply(dend, outerfuncoutermarkers)
  # attr(dend,  "nodePar")$memberlabels
  # attr(dend,  "nodePar")$outmarkers
  # 
  # 
  # attr(dend[[1]],  "nodePar")$memberlabels
  # attr(dend[[1]],  "nodePar")$outmarkers
  # 
  
  dendTR<-data.tree::as.Node(dend)
  
  # dendUnc<-dendextend::unclass_dend(dend)
  # dend <- remidpoint(dendUnc)
  # class(dend) <- "dendrogram"
  if (conditional){
    if (correction=="PC"){
      GenoCenter<-scale(SNPs,center = TRUE,scale=FALSE)
      svdGeno<-svd(GenoCenter)
      PCGeno<-GenoCenter%*%svdGeno$v
      Zall=as.matrix(dist(PCGeno[,1:nPC]))
    }
  }
  
  
  if (conditional){
    if (correction=="K"){
      GenoCenter<-scale(SNPs,center = TRUE,scale=FALSE)
      Zall=as.matrix(dist(GenoCenter))
    }
  }
  
  for (trait in 1:ncol(Pheno)){
    print(paste("Trait", trait))
    TraitData<-Pheno[,trait]
    isnotnaTrait<-!is.na(TraitData)
    x=dist(Pheno[isnotnaTrait,trait])
    
    if (conditional){
      if (correction=="PC"){
        z<-as.dist(Zall[isnotnaTrait,isnotnaTrait])
      }
      if (correction=="K"){
        z<-as.dist(Zall[isnotnaTrait,isnotnaTrait])
      }
    }
    PvalNode <- function(Node){
      pvalParent<-Node$parent$pval[trait]
      if (is.null(pvalParent)){
        pvalParent<-0
      }
      if (pvalParent<=alpha){
        print("Testing")
        inmarkers=Node$nodePar$memberlabels
        outmarkers<-Node$nodePar$outmarkers
        y=dist(SNPs[isnotnaTrait,inmarkers])
        if (length(outmarkers)>0){
          if (conditional){
            if ((correction!="PC" & correction!="K")){
              z=SNPs[isnotnaTrait,outmarkers]
            }
            #pdout<- dcorTpart.test(x,y,z,2)$p.value
           # print(pdout)
            pdout<-partialdistcov(x,y,z,1)$pdcortest
          } else {
            #pdout<- dcorT.test(x,y)$p.value
           # print(pdout)
            
            pdout<-distcov(x,y)$dcortest
          }
        } else {
         #pdout<- dcorT.test(x,y)$p.value
         #print(pdout)
         
         pdout<-distcov(x,y)$dcortest
        }
      } else {
        pdout<-1
        
        #pdout<-pvalParent
      }
      if (is.null(Node$pval)){ 
        Node$pval<-pdout
      } else {
        Node$pval<-c(Node$pval,pdout)
      }
    } 
    for (levelnode in 1:nlevels){
      print(levelnode)
      traversal <- data.tree::Traverse(dendTR, traversal = "level", filterFun = function(x){x$level==levelnode})
      data.tree::Do(traversal,PvalNode)
    }
  }
  return(dendTR)
}
