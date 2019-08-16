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

treesize.array<-c(500)#,500,800)
simtree<-200
#tol.vals.array<-rev(10^(seq(from=-16, to = -10,by=0.1)))
age.array<-10^(seq(from=12, to=16, length.out=20))
#log10kappa.array<-array(0,c(length(tol.vals.array),simtree))
log10kappa.array<-array(0,c(length(age.array),simtree))
#sol.array<-array(NA,c(length(tol.vals.array),simtree))
sol.array<-array(NA,c(length(age.array),simtree))

for(treesizeIndex in 1:length(treesize.array)){
  treesize<-treesize.array[treesizeIndex]
#  for(tolIndex in 1:(length(tol.vals.array)-1)){
#    tolmin<-tol.vals.array[tolIndex+1]
#    tolmax<-tol.vals.array[tolIndex]
  for(ageIndex in 1:length(age.array)){
    print(ageIndex)
    age<-age.array[ageIndex]
    for(reptreeIndex in 1:simtree){
#      brith_lambda=runif(1,0.01,0.1)
#      death_mu=runif(1,0,brith_lambda)
#      tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
#      tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=0.5,mu=0.1,complete = FALSE, stochsampling = TRUE)[[1]] 
      tree<-sim.bd.taxa.age(n=treesize,numbsim=1,lambda=0.5,mu=0.1,age=age,mrca=TRUE)[[1]]
      C<-vcv(tree)
      #shrink.tol<- runif(1, min=tolmin, max=tolmax)
      #C <-  tinytipvcv(C=vcv(tree),shrink=shrink.tol) 
#      log10kappa.array[tolIndex,reptreeIndex]<-log10(kappa(C))
      log10kappa.array[ageIndex,reptreeIndex]<-log10(kappa(C))
      #try(sol.array[tolIndex,reptreeIndex]<-solve(C)[1])
      try(sol.array[ageIndex,reptreeIndex]<-solve(C)[1])
      
    }
  } 

#  length(tol.vals.array)
  unsoltree<-t(apply(sol.array,1,is.na))
  unsolprop<-apply(unsoltree,1,sum)/simtree
#  names(unsolprop)<-log10(tol.vals.array)
  #print(unsolprop)
  usetable<-rbind(round(apply(log10kappa.array,1,median),2),unsolprop)
  rownames(usetable)<-c("log10 Kappa", "unsolprop")
  usetable<-usetable[,-dim(usetable)[2]]
  assign(paste("tree",treesize,sep=""),usetable)
  usetable<-NULL
}

#save.image("upperboundkappa100to800ageway.rda")
save.image("upperboundkappa500only.rda")




rm(list=ls())
load("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/upperboundkappa500.rda")
load("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/upperboundkappa00.rda")
load("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/upperboundkappa800.rda")
tree50
tree100
tree500
tree800

par(mfrow=c(2,2))
plot(tree50[1,],tree50[2,])
plot(tree100[1,],tree100[2,])
plot(tree500[1,],tree500[2,])
plot(tree800[1,],tree800[2,])
colnames(tree50)<-colnames(tree100)<-colnames(tree500)<-colnames(tree800)<-NULL
o50<-order(tree50[1,])
tree50[1,]
tree50<-tree50[,o50]
tree50[1,]
tree50[2,]
o100<-order(tree100[1,])
tree100<-tree100[,o100]
o500<-order(tree500[1,])
tree500<-tree500[,o500]
o800<-order(tree800[1,])
tree800<-tree800[,o800]


merge50<-array(NA,c(2,10))
merge100<-array(NA,c(2,10))
merge500<-array(NA,c(2,10))
merge800<-array(NA,c(2,10))

for(i in 1:10){
  merge50[1,i]<-mean(tree50[1,((i-1)*4+1):(i*4)])
  merge50[2,i]<-mean(tree50[2,((i-1)*4+1):(i*4)])
  merge100[1,i]<-mean(tree100[1,((i-1)*4+1):(i*4)])
  merge100[2,i]<-mean(tree100[2,((i-1)*4+1):(i*4)])
  merge500[1,i]<-mean(tree500[1,((i-1)*4+1):(i*4)])
  merge500[2,i]<-mean(tree500[2,((i-1)*4+1):(i*4)])
  merge800[1,i]<-mean(tree800[1,((i-1)*4+1):(i*4)])
  merge800[2,i]<-mean(tree800[2,((i-1)*4+1):(i*4)])
  }

library(ggplot2)
df<-data.frame( log10kappa=c(merge50[1,],merge100[1,],merge500[1,],merge800[1,]),
                Proportion=c(1-merge50[2,],1-merge100[2,],1-merge500[2,],1-merge800[2,]), 
                Taxa=c(rep(50,dim(merge50)[2]),rep(100,dim(merge50)[2]),rep(500,dim(merge50)[2]),rep(800,dim(merge50)[2]))
)

df$Taxa<-as.factor(df$Taxa)
p.plot<-ggplot(df,aes(x=log10kappa,y=Proportion,col=Taxa,group=Taxa))+ geom_point() + geom_line(size=2)
p.plot<- p.plot + ggtitle(label = "Solvable C Proportion") + labs(x=expression(paste(log[10], " ", kappa))) #+ xlim(14,20) 
p.plot<-p.plot+theme(
  plot.title = element_text(color = "black", size = 20, face = "bold"))+ theme(plot.title = element_text(hjust = 0.5))

p.plot


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ggtree")
library(ape)
library(ggtree)
treesize<-6
tree<-sim.bd.taxa.age(n=treesize,numbsim=1,lambda=0.8, mu=0.1,age=100,mrca=TRUE)[[1]]
tip.label<-c("  A","  B","  C","  D","  E","  F")
tree$tip.label<-tip.label
plot(tree,edge.width=3)


nodelabels("",1,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",2,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",3,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",4,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",5,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",6,frame="c",bg="black",adj=0,cex=0.5)


axisPhylo(1,las=1,backward=FALSE)
#edgelabels(round(tree$edge.length,2),bg="black",adj=c(0.5,0.5),col="white",font=2)


