rm(list=ls())
#simulate a tree of 10 taxa 
#install.packages("TreeSim")
library(xtable)
library(TreeSim)
size<-6
trees <-  sim.bd.taxa.age(n=size, numbsim=1, lambda=0.1, mu=0.1, frac = 1, age=100, mrca = TRUE)#lambda:speciation rate, mu:extinction rate, frac: each tip is included into the final tree with prob frac

phy<-trees[[1]]
phy$tip.label<-c("  A","  B","  C","  D","  E","  F")
plot(phy,edge.width=5,xex=1.5)

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
edgelabels(round(phy$edge.length,2),bg="black",adj=c(0.5,0.5),col="white",font=2)


xtable(round(vcv(phy),2))