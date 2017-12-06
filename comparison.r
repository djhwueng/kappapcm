#idea: simulate many trees under a given rate, estimate BM rate under a variety of transforms(removing tips, shrinkage method, eigenvalue method) for the matrix. Look at the RMSE for BM rate under each of these approaches. Presumably one will win. And we can code that in BMhyb package or, even make a separate mvn package that calculates likelihood given a vcv, mu vector, and states, and automatically does a transforms.


#test rate RMSE
#test rate RMSE

#test rate RMSE
rm(list=ls())
setwd("/Users/djhwueng/Dropbox/Collab")
library(phyclust)
library(geiger)
library(TreeSim)
library(ape)
library(Matrix)
library(phytools)
library(phangorn)
library(geiger)
library(MASS)

getShortest<-function(phy){
  tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  tip.lengths
  return(which.min(tip.lengths))
}

DropShortestTipTrait<-function(phy,traitset){
  tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  drop.tip.index<-phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]]
  drop.tip.names<-names(tip.lengths[drop.tip.index])
  phy<-drop.tip(phy, drop.tip.index)
  traitset<-traitset[- which( rownames(traitset)== drop.tip.names)  ,]
  return(list(phy=phy,traitset=traitset))
}


DropRandom<-function(phy){
  phy<-drop.tip(phy,sample.int(Ntip(phy),1))
  return(phy)
}

omega<-function(phy){
  tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  tree.height<-get.rooted.tree.height(phy)
  return(min(tip.lengths)/tree.height)
}

illcond.tree<-function(phy,eta){
  return(omega(phy)<eta)
}

#Theiler's shrinkage estimator (2012)
TH.cov.shrink<-function(S,n){
  T<-diag(diag(S))
  alpha<-10^(optimize(avenegLOOL,c(-8,0),S,T,n)$minimum)
  sig.pk<-(1-alpha)*S + alpha*T
  sig.pk
}

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

shrink.tree<-function(phy){
  phy.vcv<-vcv(phy)
  vcv.shrunk<-TH.cov.shrink(phy.vcv,dim(phy.vcv)[1])
  shrunk.tree<-upgma(2*(max(vcv.shrunk)-vcv.shrunk))
  return(shrunk.tree)
}

GetAveSigsq<-function(one.tree.result){
  return(one.tree.result$opt$sigsq)
}

rmse<-function(obs,true=true){
  sqrt(mean((obs-true)^2))
}


sims<-50 #both replicates of  tree and trait
sigsq<-1
treesize<-100
eta<-1e-2#checking tree condition
#trees<-sim.bd.taxa.age(treesize,sims,lambda=0.4,mu=0.1,frac=0.5,age=log(100)/0.3)
trees <-  sim.bd.taxa.age(n=treesize, numbsim=sims, lambda=1, mu=0.5, frac = 0.5, age=100, mrca = TRUE)

results<-list(raw=NULL,drop=NULL,shrink=NULL)
num.illtree<-0
for(treeIndex in 1:sims){
  phy<-trees[[treeIndex]]
  if(illcond.tree(phy,eta)){#we have ill cond tree
    num.illtree<-num.illtree+1
    phy.raw<-phy
    traitset<-sim.char(phy.raw,par=sigsq,nsim=sims,model="BM",root=1)
    traitset<-drop(traitset)
    rownames(traitset)<-phy$tip.label

    fitC.raw<-apply(traitset, 2, fitContinuous, phy=phy.raw, model="BM", control = list(niter=10))


    dropping<- DropShortestTipTrait(phy.raw,traitset)
    phy.drop<-dropping$phy
    traitset.drop<-dropping$traitset

    fitC.drop<-apply(traitset.drop, 2, fitContinuous, phy=phy.drop, model="BM", control = list(niter=10))

    phy.shrink<-shrink.tree(phy)
    fitC.shrink<-apply(traitset, 2, fitContinuous, phy=phy.shrink, model="BM", control = list(niter=10))

    temp.raw<-matrix(sapply(fitC.raw,GetAveSigsq),ncol=1)
    results$raw<-cbind(results$raw,temp.raw)

    temp.drop<-matrix(sapply(fitC.drop,GetAveSigsq),ncol=1)
    results$drop<-cbind(results$drop,temp.drop)

    temp.shrink<-matrix(sapply(fitC.shrink,GetAveSigsq),ncol=1)
    results$shrink<-cbind(results$shrink,temp.shrink)

  }
}
print(num.illtree)
print(sapply(results,rmse,true=sigsq))

save.image("comparison.RData")


#plot rmse vs kappa to check whether linear relationship exists.
rm(list=ls())
setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara")
load("comparison.RData")
kappaget<-function(phy){
  kappa(vcv(phy))
}

par(mfrow=c(1,3))
plot( log(sapply(trees,kappaget)),apply(results$raw,2,rmse,true=sigsq),xlab="log(kappa)",ylab="rmse for sigms.sq", main="Raw Tree")
plot( log(sapply(trees,kappaget)),apply(results$drop,2,rmse,true=sigsq),xlab="log(kappa)",ylab="rmse for sigms.sq", main="Drop the Shortest tip Tree")
plot( log(sapply(trees,kappaget)),apply(results$shrink,2,rmse,true=sigsq), ,xlab="log(kappa)",ylab="rmse for sigms.sq", main="Shrunk Tree")



summary(lm(apply(results$raw,2,rmse,true=sigsq) ~ log(sapply(trees,kappaget))))

summary(lm(apply(results$drop,2,rmse,true=sigsq) ~ log(sapply(trees,kappaget))))
  
summary(lm(apply(results$shrink,2,rmse,true=sigsq) ~ log(sapply(trees,kappaget))))

