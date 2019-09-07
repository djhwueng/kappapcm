setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaShrinkRun3/")
rm(list=ls())
#install.packages("tidyr")
#install.packages("pracma")
#install.packages("SEAsic")

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
#library(SEAsic)



TH.cov.shrink.delta<-function(S,n,delta){
  T<-diag(diag(S))
  #  delta<-10^(optimize(avenegLOOL,c(-10,0),S,T,n)$minimum)
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


Shrink.C<-function(C,delta){
  phy.vcv<-C
  vcv.shrunk<-TH.cov.shrink.delta(phy.vcv,dim(phy.vcv)[1],delta)
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



RMSD_sametrait<-function(Yset=Yset,adjC=adjC,droptipname=NULL,method=method){
  sims<-dim(Yset)[2]
  param.array<-array(0,c(sims,2))
  for (i in 1:sims) {
    y<-Yset[,i]
    names(y)<-rownames(adjC)
    if(!is.null(droptipname)){
      y<-y[-which((names(y)==droptipname)) ]
    }
    param.array[i,]<- MLEfcn(adjC=adjC,Y=y,method=method)
  }
  colnames(param.array)<-c("theta","sigmasq")
  rmsd_mu<-  sqrt(mean((param.array[,1]-mu) ,na.rm=TRUE)^2)    #sqrt(mean((param.array[,1]-mu)^2))/mean(param.array[,1],na.rm=TRUE)
  rmsd_sigmasq <- sqrt(mean(param.array[,2]-sigmasq, na.rm=TRUE)^2)#sqrt(mean((param.array[,2]-sigmasq)^2))/mean(param.array[,2],na.rm=TRUE)
  return( c(rmsd_mu,rmsd_sigmasq))
}


############# Brownian model
MLEfcn<-function(adjC=adjC,Y=Y,method=method){
  one=matrix(1,ncol=1,nrow=length(Y))
  if (method=="solve"){solveM<-solve(adjC,tol=1e-100)} # HERE IS THE KEY TO AFFECT KAPPA AGAIN
  if (method=="svd"){solveM<-svd_decom(adjC)}
  if (method=="ginv"){solveM<-ginv(adjC,tol=1e-100)}
  mu<-c((t(one)%*%solveM%*%Y)/(t(one)%*%solveM%*%one))
  Y<-Y-mu
  sigma_sq<-c((t(Y)%*%solveM%*%Y)/length(Y))
  return(c(mu,abs(sigma_sq)))
}

# Idea simlates trees

simsTree<-100# that means repeat 100 times
treesize<-100 # 100, 500, 800
simtrait<-50 # 50 traits
Yset<-array(0,c(treesize,simtrait))
mu<-10 # or mu <- 10
sigmasq<-5# sigmasq<- 5

savename<-paste("unsolshrinksize",treesize,"mu",mu,"sigmasq",sigmasq,".rda",sep="")
savename

delta.array<- c(0.001,0.005,0.01,1)
#solve_raw_RMSD_combined<-NULL
for(deltaIndex in 1:length(delta.array)){
  print(deltaIndex)
  delta<-delta.array[deltaIndex]
  #solve_raw_RMSD<-NULL
  
  unsolve_shrink_RMSD<-NULL
  for(treeIndex in 1:simsTree){
    print(treeIndex)
    brith_lambda=runif(1,0.01,0.1)
    death_mu=runif(1,0,brith_lambda)
    tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
    C<-vcv(tree)
    C<-C/max(C)
    kappaC<-kappa(C)
    shrink.tol<- runif(1, min=1e-15, max=5e-13)
    C <-  tinytipvcv(C=vcv(tree),shrink=shrink.tol) # Bad Tree we use a very bad tree to simulate data, as expected parameter estimate will be bad. so adjust tree will lead a better estimate. need to determine the bound for the tree.
    #print(kappa(C))  
    one<-array(1,c(dim(C)[1],1))
    for (simtraitIndex in 1:simtrait) {
      Yset[,simtraitIndex]<- mvrnorm(n = 1, mu*one, Sigma=sigmasq*C, tol = 1e-10, empirical = FALSE, EISPACK = FALSE)
    }
    #try(solve_raw_RMSD<-rbind(solve_raw_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
    C<-Shrink.C(C,delta)
    try(unsolve_shrink_RMSD<-rbind(unsolve_shrink_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  }
  # solve_raw_RMSD_combined<-rbind(solve_raw_RMSD_combined,solve_raw_RMSD)
  assign(paste("unsolve_shrink_RMSD_delta",delta,sep=""),unsolve_shrink_RMSD)
  save.image(savename)
}



