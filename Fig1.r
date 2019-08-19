rm(list=ls())
#simulate a tree of 10 taxa 
#install.packages("TreeSim")
library(xtable)
library(TreeSim)


tinytipvcv<-function(C=C,shrink=shrink){
  #C<-C/max(C)
  uniC<-unique(c(C))
  secondlargest<-  sort(uniC)[(length(uniC)-1)]
  diag(C)<-(max(C)-secondlargest)*shrink+secondlargest
  return(C/max(C))
  #return(C)
}

set.seed(58)
size<-6
tree <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=0.6, mu=0.1, frac = 1, age=1, mrca = TRUE)[[1]]#lambda:speciation rate, mu:extinction rate, frac: each tip is included into the final tree with prob frac
#tree<-sim.bd.taxa(n=size,numbsim=1,lambda=0.7,mu=0.5,complete = FALSE, stochsampling = TRUE)[[1]] 
tree$tip.label<-c("  A","  B","  C","  D","  E","  F")
plot(tree,edge.width=5)

nodelabels("",1,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",2,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",3,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",4,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",5,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",6,frame="c",bg="black",adj=0,cex=0.5)

nodelabels("",7,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("G",7,frame="c",bg="black",adj=-1.5,cex=1)
nodelabels("",8,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("H",8,frame="c",bg="black",adj=-1.5,cex=1)
nodelabels("",9,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("I",9,frame="c",bg="black",adj=-2.5,cex=1)
nodelabels("",10,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("J",10,frame="c",bg="black",adj=-1.5,cex=1)
nodelabels("",11,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("K",11,frame="c",bg="black",adj=-1.5,cex=1)

axisPhylo(1,las=1,backward=FALSE)
edgelabels(round(tree$edge.length,2),bg="black",adj=c(0.5,0.5),col="white",font=2)


xtable(round(vcv(tree),2))




par(mfrow = c(1,2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
# tiny tree
library(phangorn)
C<-vcv(tree)
C <- tinytipvcv(C=vcv(tree),shrink=10^(-3))
tree1<-upgma(2*(max(C)-C))
tip.label<-c("  A","  B","  C","  D","  E","  F")
tree1$tip.label<-tip.label
plot(tree1,edge.width=5)


nodelabels("",1,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",2,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",3,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",4,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",5,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",6,frame="c",bg="black",adj=0,cex=0.5)


axisPhylo(1,las=1,backward=FALSE)


tree2<-sim.bd.taxa.age(n=size,numbsim=1,lambda=0.6,mu=0.1,age=9*10^2,mrca=TRUE)[[1]] 
C <- vcv(tree2)
C<- C/max(C)
tree2<-upgma(2*(max(C)-C))
tip.label<-c("  A","  B","  C","  D","  E","  F")
tree2$tip.label<-tip.label
plot(tree2,edge.width=5)
nodelabels("",1,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",2,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",3,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",4,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",5,frame="c",bg="black",adj=0,cex=0.5)
nodelabels("",6,frame="c",bg="black",adj=0,cex=0.5)


axisPhylo(1,las=1,backward=FALSE)

print(c(kappa(vcv(tree1)),kappa(vcv(tree2))))

