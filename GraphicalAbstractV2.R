rm(list=ls())
#simulate a tree of 10 taxa 
#install.packages("TreeSim")
library(phyclust)
library(geiger)
library(TreeSim)
library(ape)
library(Matrix)
library(phytools)
library(phangorn)
library(MASS)


#Theiler 2012 optimization to find alpha
avenegLOOL<-function(logalpha,S,T,n){
  alpha<-10^logalpha; beta<-(1-alpha)/(n-1)
  G<-beta*n*S +alpha*T
  elu<-expand(lu(G))
  L<-elu$P%*%elu$L; U<-elu$U
  logGdet<-sum(log(abs(diag(U))))
  ro<-sum(diag(solve(G)%*%S))
  if(beta*ro>1){
    val<-Inf
  }else{
    val<-log(1-beta*ro)+ro/(1-beta*ro) + logGdet
  }
  return(val)
}


#Theiler's shrinkage estimator (2012)
TH.cov.shrink<-function(S,n){
  T<-diag(diag(S))
  alpha<-10^(optimize(avenegLOOL,c(-8,0),S,T,n)$minimum)
  sig.pk<-(1-alpha)*S + alpha*T
  sig.pk
}

shrink.tree<-function(phy){
  phy.vcv<-vcv(phy)
  vcv.shrunk<-TH.cov.shrink(phy.vcv,dim(phy.vcv)[1])
  shrunk.tree<-upgma(2*(max(vcv.shrunk)-vcv.shrunk))
  return(shrunk.tree)
}

DropShortestTip<-function(phy){
  tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  drop.tip.index<-phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]]
  drop.tip.names<-names(tip.lengths[drop.tip.index])
  phy<-drop.tip(phy, drop.tip.index)
  return(phy)
}


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
 return(round(kappa(vcv(phy)),1))
}

set.seed(1623)
size<-6
trees <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=0.1, mu=0.1, frac = 1, age=100, mrca = TRUE)#lambda:speciation rate, mu:extinction rate, frac: each tip is included into the final tree with prob frac


par(mfrow=c(2,2))

op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,2,2) + 0.1)


phy<-trees[[1]]
phy$tip.label<-letters[1:size]
plot(phy,main=paste("Raw Tree: k = ", kappaget(phy),sep=""),edge.width=4,show.tip.label=FALSE,direction="upwards",cex.main=2)
phy.d<-DropShortestTip(phy)
plot(phy.d, main=paste("Pruned Tree: k = ",kappaget(phy.d),sep=""),edge.width=4,show.tip.label=FALSE,direction="upwards",cex.main=2)
phy.s<-shrink.tree(phy)
plot(phy.s,main=paste("Shrunk Tree: k = ",kappaget(phy.s),sep=""),edge.width=4,show.tip.label=FALSE,direction="upwards",cex.main=2)
phy.p<-PerturbTree_beta(phy)
plot(phy.p,main=paste("Stretched Tree: k = ",kappaget(phy.p),sep=""),edge.width=4,show.tip.label=FALSE,direction="upwards",cex.main=2)
par(op)
