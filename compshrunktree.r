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
library(ggplot2)


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
  #alpha<-10^(optimize(avenegLOOL,c(-8,0),S,T,n)$minimum)
  alpha<-0.5
  sig.pk<-(1-alpha)*S + alpha*T
  sig.pk
}

shrink.tree<-function(phy){
  phy.vcv<-vcv(phy)
  vcv.shrunk<-TH.cov.shrink(phy.vcv,dim(phy.vcv)[1])
  shrunk.tree<-upgma(2*(max(vcv.shrunk)-vcv.shrunk))
  return(shrunk.tree)
}

set.seed(5812)
size<-100
#phy <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=0.1, mu=0.1, frac = 1, age=10, mrca = TRUE)[[1]]#lambda:speciation rate, mu:extinction rate, frac: each tip is included into the final tree with prob frac
phy <-  sim.bd.taxa(n=size, numbsim=1, lambda=0.5, mu=0.2, complete=FALSE)[[1]]#lambda:speciation rate, mu:extinction rate, frac: each tip is included into the final tree with prob frac
#phy <-  rcoal(size)
phy.s<-shrink.tree(phy)
phy$edge
phy.s$edge

par(mfrow = c(1, 2),     # 2x2 layout
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA)            # allow content to protrude into outer margin (and beyond)

phy.dat<-cbind(phy$edge,phy$edge.length)
phy.s.dat<-cbind(phy.s$edge,phy.s$edge.length)
phy.s2.dat<-array(NA,dim(phy.s.dat))
#phy.s2.tip.label<-phy.s$tip.label
for(i in 1:dim(phy.dat)[1]){
  #print(which(phy.dat[,2] %in% phy.s.dat[i,2]))
  print(phy.dat[i,2])
  matchindex<- which(phy.s.dat[,2] %in% phy.dat[i,2])
  phy.s2.dat[i,]<-phy.s.dat[matchindex,] 
  phy.dat
  phy.s2.dat
}

phy.dat
phy.s2.dat
phy.s$edge.length<-phy.s2.dat[,3]
phy.s$edge<-phy.s2.dat[,1:2]
phy.s$edge
phy$edge
plot(phy,show.tip.label = FALSE,edge.width = 2)
plot(phy.s,show.tip.label = FALSE,edge.width = 2)






