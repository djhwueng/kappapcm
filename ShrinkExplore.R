rm(list=ls())
setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationSummary")
rm(list=ls())

library(tidyr)
library(ape)
library(MASS)
library(TreeSim)
library(Matrix)
library(phytools)
library(phangorn)
library(dplyr)
library(corpcor)
library(pracma)

############## shrink method
avenegLOOL<-function(logdelta,C,T,n){
  delta<-10^logdelta
  beta <-(1-delta)/(n-1)
  S_delta<-n*beta*C+delta*T #formula (3) in
  #  elu<-expand(lu(S_delta))
  elu<-lu(S_delta)
  #elu<-lu(S_delta)
  #L<-elu$P%*%elu$L
  L<-elu$L
  U<-elu$U
  logSdet<-sum(log(abs(diag(U))))
  ro<-sum(diag(solve(S_delta)%*%C))
  if(beta*ro>1){
    val<-Inf
  }else{
    val<-log(1-beta*ro)+ro/(1-beta*ro)+logSdet
  }
  return(val)
}


TH.cov.shrink<-function(S,n){
  T<-diag(diag(S))
  delta<-10^(optimize(avenegLOOL,c(-10,0),S,T,n)$minimum)
  beta <-(1-delta)/(n-1)
  C<-n*beta*S + delta*T
  return(C)
}

Shrink.tree<-function(phy){
  phy.vcv<-vcv(phy)
  vcv.shrunk<-TH.cov.shrink(phy.vcv,dim(phy.vcv)[1])
  shrunk.tree<-upgma(2*(max(vcv.shrunk)-vcv.shrunk))
  return(shrunk.tree)
}


Shrink.C<-function(C){
  phy.vcv<-C
  vcv.shrunk<-TH.cov.shrink(phy.vcv,dim(phy.vcv)[1])
  return(vcv.shrunk/max(vcv.shrunk))
}


# NEED to check tinytipvcv MAKE THIS MATRI POSTIVE DEFINTE
tinytipvcv<-function(C=C,shrink=shrink){
  #C<-C/max(C)
  uniC<-unique(c(C))
  secondlargest<-  sort(uniC)[(length(uniC)-1)]
  diag(C)<-(max(C)-secondlargest)*shrink+secondlargest
  return(C/max(C))
  #return(C)
}


simsTree<-500
#treesize.array<-c(20,100,150,500, 800) 
treesize.array<-c(100) 
beta.delta.array <- array(0,c(length(treesize.array),simsTree,2))

for(treesizeIndex in 1: length(treesize.array)){
  treesize<- treesize.array[treesizeIndex]
  for(treeIndex in 1:simsTree){
    print(treeIndex)
    brith_lambda=runif(1,0.01,0.1)
    death_mu=runif(1,0,brith_lambda)
    tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
    shrink.tol<- runif(1, min=1e-15, max=5e-13)
  #  C <- tinytipvcv(C=vcv(tree),shrink=shrink.tol) # Bad Tree we use a very bad tree to simulate data, as expected parameter estimate will be bad. so adjust tree will lead a better estimate. need to determine the bound for the tree.
    C<-vcv(tree)
    T <- diag(diag(C))
    n<-dim(C)[1]
    delta<-10^(optimize(avenegLOOL,c(-10,0),C,T,n)$minimum)
    beta <-(1-delta)/(n-1)
    print(c(n*beta,delta))  
    beta.delta.array[treesizeIndex,treeIndex,]<-c(n*beta,delta)
  }
}


for(Index in 1:length(treesize.array)){
  print(paste(round(apply(beta.delta.array[Index,,],2,mean),2)," treesize =", treesize.array[Index],sep=""))
  }
  system("pwd")  
#save.image("ShrinkExplore.rda")
#load("ShrinkExplore.rda")



