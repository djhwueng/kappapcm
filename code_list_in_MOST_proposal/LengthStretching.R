rm(list=ls())
#----------------------------------------------
SegM<-function(M){#find branch length segments
  segment<-NULL
  sortM<-sort(M)
  while(sum(sortM!=0)){ #from tip to bottom 
    m<-max(sortM)        
    sortM[which(sortM==m)]<-0
    sm<-max(sortM)
    segment<-cbind(segment,m-sm)# t_d,...,t_2,t_1
    }
  return(rev(segment)) #t_1,t_2,...,t_d from root to tips.
  }

PtbLen<-function(sgmt){ #varying time segments
  k<-ceiling(1/min(sgmt)) 
  applen<-array(0,c(length(sgmt)))
  for (i in 1:length(sgmt) ){
    applen[i]<-rgamma(1,shape=k*sgmt[i],rate=1) 
    }
  return(applen/sum(applen)) #Dirichlet samples
  }
 

PtbLen_beta<-function(sgmt){ #varying time segments under beta
  new_sgmt<-array(0,length(sgmt))
  new_sgmt[1]<-rbeta(n=1,shape1 = sgmt[1],shape2 = sum(sgmt[-1]))
  for(Index in 2:(length(sgmt)-1)){
    s_Index<-rbeta(n=1,shape1 = sgmt[Index],shape2 = sum(sgmt[-(1:Index)]))
    new_sgmt[Index]<-(1-sum(sgmt[1:(Index-1)]))*s_Index
  }
  new_sgmt[length(sgmt)] <- 1 - sum(new_sgmt[-length(new_sgmt)]) 
  return(new_sgmt)
}




AccSum<-function(sgmt){ #accumulate segments
  accsum<-array(0,c(length(sgmt))) 
  tmp<-0
  for(i in 1:length(sgmt)){
    tmp<-tmp+sgmt[i]   
    accsum[i]<-tmp
    }
  return(accsum)
  } #g1=t1,g2=t1+t2,...,gd=t1+t2+...+td
  

VOR<-function(accsum, mtx=mtx){ #obtain the matrix
  accsum<-rev(accsum) #gd, ... ,g2 ,g1
  outmtx<-array(0,c(dim(mtx)))
  M<-mtx
  count<-0
  while(sum(M)!=0){
    count<-count+1
    ptb<-max(M)
    positions<-which(mtx==ptb)
    outmtx[positions]<-accsum[count]   
    M[positions]<-0  
    }
  return(outmtx)
  }  
  
PerturbTree<-function(phy){ #get the new tree
  mtx<-vcv(phy)
  height<-max(mtx)
  mtx<-mtx/height #scale tree to make the Dirichlet distribution work
  sgmt<- SegM(mtx) #get segments
  ptblen<-PtbLen(sgmt) #perturb segments
  accsum<-AccSum(ptblen)
  M<-VOR(accsum,mtx=mtx)
  rownames(M)<-rownames(mtx)
  colnames(M)<-colnames(mtx)
  M<-height*M # put the tree height back for the perturb tree
  phy<-upgma(2*(max(M)-M))
  return(phy)
  }


PerturbTree_beta<-function(phy){ #get the new tree
  mtx<-vcv(phy)
  height<-max(mtx)
  mtx<-mtx/height #scale tree to make the Dirichlet distribution work
  sgmt<- SegM(mtx) #get segments
  ptblen<-PtbLen_beta(sgmt) #perturb segments
  accsum<-AccSum(ptblen)
  M<-VOR(accsum,mtx=mtx)
  rownames(M)<-rownames(mtx)
  colnames(M)<-colnames(mtx)
  M<-height*M # put the tree height back for the perturb tree
  phy<-upgma(2*(max(M)-M))
  return(phy)
}

kappaget<-function(phy){
  kappa(vcv(phy))
  }

#library(phyclust)
#library(geiger)
library(TreeSim)
library(ape)
#library(Matrix)
#library(phytools)
library(phangorn)
#library(geiger)
#library(MASS)

 
treesize<-100
numbsim<-1 #simulate 1 tree, treat it as the true tree.
trees<-sim.bd.taxa.age(n=treesize,numbsim=numbsim,lambda=1,mu=0.5, frac=0.5, age= (log(treesize)/(1-0.5)))

pKappaset<-replicate(n=100,kappaget(PerturbTree(trees[[1]]))) # perturb the tree by stretching the lengths. obtain 100 kappa.  
pKappaset_beta<-replicate(n=100,kappaget(PerturbTree_beta(trees[[1]]))) # perturb the tree by stretching the lengths. obtain 100 kappa.  

par(mfrow=c(1,2))
plot(log(pKappaset),main="Length Stretch: Gamma",ylab="log kappa", xlab="replicates index")
abline(h=log(kappa(vcv(trees[[1]]))),col="red")
plot(log(pKappaset_beta),main="Length Stretcg: Beta",ylab="log kappa", xlab="replicates index")
abline(h=log(kappa(vcv(trees[[1]]))),col="red")


sort()

par(mfrow=c(1,3))
plot(trees[[1]],main="raw tree")
plot(PerturbTree(trees[[1]]),main="new tree:gamma")
plot(PerturbTree_beta(trees[[1]]),main="new tree:beta")
