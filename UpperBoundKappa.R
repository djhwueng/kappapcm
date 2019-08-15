rm(list=ls())
setwd("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/")
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
# NEED to check tinytipvcv MAKE THIS MATRI POSTIVE DEFINTE
tinytipvcv<-function(C=C,shrink=shrink){
  #C<-C/max(C)
  uniC<-unique(c(C))
  secondlargest<-  sort(uniC)[(length(uniC)-1)]
  diag(C)<-(max(C)-secondlargest)*shrink+secondlargest
  return(C/max(C))
  #return(C)
}

treesize.array<-c(50,100,500,800)
simtree<-500
tol.vals.array<-rev(10^(seq(from=-16, to = -10,by=0.1)))
log10kappa.array<-array(0,c(length(tol.vals.array),simtree))
sol.array<-array(NA,c(length(tol.vals.array),simtree))

for(treesizeIndex in 1:length(treesize.array)){
  treesize<-treesize.array[treesizeIndex]
  for(tolIndex in 1:(length(tol.vals.array)-1)){
    tolmin<-tol.vals.array[tolIndex+1]
    tolmax<-tol.vals.array[tolIndex]
    for(reptreeIndex in 1:simtree){
      brith_lambda=runif(1,0.01,0.1)
      death_mu=runif(1,0,brith_lambda)
      tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
      shrink.tol<- runif(1, min=tolmin, max=tolmax)
      C <-  tinytipvcv(C=vcv(tree),shrink=shrink.tol) 
      log10kappa.array[tolIndex,reptreeIndex]<-log10(kappa(C))
      try(sol.array[tolIndex,reptreeIndex]<-solve(C)[1])
    }
  } 
  
  length(tol.vals.array)
  unsoltree<-t(apply(sol.array,1,is.na))
  unsolprop<-apply(unsoltree,1,sum)/simtree
  names(unsolprop)<-log10(tol.vals.array)
  #print(unsolprop)
  usetable<-rbind(round(apply(log10kappa.array,1,median),2),unsolprop)
  rownames(usetable)<-c("log10 Kappa", "unsolprop")
  usetable<-usetable[,-dim(usetable)[2]]
  assign(paste("tree",treesize,sep=""),usetable)
  usetable<-NULL
}

#save.image("upperboundkappa.rda")


rm(list=ls())
#load("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/upperboundkappa500.rda")
#load("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/upperboundkappa800.rda")

load("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/upperboundkappa.rda")
par(mfrow=c(2,2))
plot(tree50[1,],tree50[2,])
plot(tree100[1,],tree100[2,])
plot(tree500[1,],tree500[2,])
plot(tree800[1,],tree800[2,])



library(ggplot2)
df<-data.frame( log10kappa=c(tree50[1,],tree100[1,],tree500[1,],tree800[1,]),
                Proportion=c(1-tree50[2,],1-tree100[2,],1-tree500[2,],1-tree800[2,]), 
                Taxa=c(rep(50,dim(tree50)[2]),rep(100,dim(tree50)[2]),rep(500,dim(tree50)[2]),rep(800,dim(tree50)[2]))
               )
df$Taxa<-as.factor(df$Taxa)
p.plot<-ggplot(df,aes(x=log10kappa,y=Proportion,col=Taxa,group=Taxa)) + geom_point() + geom_line(size=1.5)
p.plot<- p.plot + ggtitle(label = "Solvable C") + labs(x="log10 kappa")
p.plot<-p.plot+theme(
  plot.title = element_text(color = "black", size = 20, face = "bold"))+ theme(plot.title = element_text(hjust = 0.5))

p.plot

