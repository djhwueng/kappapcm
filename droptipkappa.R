setwd("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/")
rm(list=ls())
library(TreeSim)
library(ape)

DropShortesttip<-function(phy){
  tip.lengths <-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  drop.tip.index<-phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]]
  drop.tip.names<-names(tip.lengths[drop.tip.index])
  phy<-drop.tip(phy,drop.tip.index)
  return(list(phy=phy,drop.tip.names=drop.tip.names))
}

sims<-100
tree.log10kappa<-array(0,c(sims))
tree.r.log10kappa<-array(0,c(sims))
tree.s.log10kappa<-array(0,c(sims))

for(i in 1:sims){
  ntax<-runif(1,min=100,max=800)
  lambda=runif(1,0.01,0.1)
  mu=runif(1,0,lambda)
  tree<-sim.bd.taxa(n=ntax,numbsim=1,lambda=lambda,mu=mu,complete = FALSE, stochsampling = TRUE)[[1]]
  
  tree.log10kappa[i]<-log10(kappa(vcv(tree)))
  tree.r<-drop.random(phy=tree,n=1)
  tree.r.log10kappa[i]<-log10(kappa(vcv(tree.r)))
  tree.s<-DropShortesttip(phy=tree)$phy
  tree.s.log10kappa[i]<-log10(kappa(vcv(tree.s)))
}
png("fig6.png")
plot(tree.log10kappa,tree.r.log10kappa, pch=0, col="blue",
     main="Drop shortest tip vs. Original", ylab= expression(paste(" Drop 1 tip: ", log[10], kappa ,   sep="  ") )
     ,xlab= expression(paste(" Original: ", log[10], kappa ,   sep="  ") )
     )
points(tree.log10kappa,tree.s.log10kappa,col="red",pch=16)
segments(tree.r.log10kappa,tree.r.log10kappa,tree.r.log10kappa,tree.s.log10kappa )
dev.off()