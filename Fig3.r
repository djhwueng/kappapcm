setwd("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/")
rm(list=ls())
#install.packages("repmis")
library(repmis)
library(ape)
library(TreeSim)
library(phytools)
kappaget<-function(phy){
  kappa(vcv(phy))
}
treesize.array<-c(10,seq(50,800,by=50))
numbsim<-100
kappa.array<-array(0,c(numbsim,length(treesize.array)))#sim.bd.taxa.age method
kappa.array1<-array(0,c(numbsim,length(treesize.array)))#rcoal method
kappa.array2<-array(0,c(numbsim,length(treesize.array)))#ptree method

minegn.obs.array<-NULL
minegn.obs1.array<-NULL
minegn.obs2.array<-NULL

trees1<-list()
trees2<-list()
for(treesizeIndex in 1:length(treesize.array)){
  print(treesize.array[treesizeIndex])
  #birth-death tree
  #trees<- sim.bd.taxa.age(n=treesize.array[treesizeIndex],numbsim=numbsim,lambda=0.4,mu=0.1,frac=.5,age=log(100)/0.3)
  tree.size<-treesize.array[treesizeIndex]
  brith_lambda=runif(1,0.01,0.1)
  death_mu=runif(1,0,brith_lambda)
  #frac=runif(1,0.1,1)
  #age=log(tree.size)/(brith_lambda-death_mu)
  #trees<-sim.bd.taxa.age(n=tree.size,numbsim=numbsim,lambda=brith_lambda,mu=death_mu,frac=frac,age=log(treesize.array[treesizeIndex])/(brith_lambda-death_mu))
  
  #trees<-sim.bd.taxa.age(n=tree.size,numbsim=numbsim,lambda=brith_lambda,mu=death_mu,frac=frac,age=1)
  
  trees<-sim.bd.taxa(n=tree.size,numbsim=numbsim,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)
  
  vcv.obs<-lapply(1:numbsim,function(x) vcv(trees[[x]])/diag(vcv(trees[[x]]))[1])
  minegn.obs<-lapply(1:numbsim,function(x)  min(eigen(vcv(trees[[x]])/diag(vcv(trees[[x]]))[1])$values))
  minegn.obs.array <- cbind(minegn.obs.array, matrix(unlist(minegn.obs),ncol=1))
  
  
  kappa.obs<-unlist(lapply(1:numbsim, function(x) kappa(vcv.obs[[x]],exact=TRUE)))
  kappa.array[,treesizeIndex]<-kappa.obs
  
  #random split tree
  for(Index in 1:numbsim){
    trees1[[Index]]<-rcoal(treesize.array[treesizeIndex])
  }
  
  vcv.obs1<-lapply(1:numbsim, function(x) vcv(trees1[[x]])/diag(vcv(trees1[[x]]))[1])
  minegn.obs1<-lapply(1:numbsim,function(x)  min(eigen(vcv(trees1[[x]])/diag(vcv(trees1[[x]]))[1])$values))
  minegn.obs1.array <- cbind(minegn.obs1.array, matrix(unlist(minegn.obs1),ncol=1))
  
  kappa.obs1<-unlist(lapply(1:numbsim, function(x) kappa(vcv.obs1[[x]],exact=TRUE)))
  kappa.array1[,treesizeIndex]<-kappa.obs1
  
  #pure birth tree
  for(Index in 1:numbsim){
    trees2[[Index]]<-pbtree(n=treesize.array[treesizeIndex])
  }
  vcv.obs2<-lapply(1:numbsim, function(x) vcv(trees2[[x]])/diag(vcv(trees2[[x]]))[1])
  minegn.obs2<-lapply(1:numbsim,function(x)  min(eigen(vcv(trees2[[x]])/diag(vcv(trees2[[x]]))[1])$values))
  minegn.obs2.array <- cbind(minegn.obs2.array, matrix(unlist(minegn.obs2),ncol=1))
  
  kappa.obs2<-unlist(lapply(1:numbsim, function(x) kappa(vcv.obs2[[x]],exact=TRUE)))
  kappa.array2[,treesizeIndex]<-kappa.obs2
}

summary(kappa.array)
summary(kappa.array1)
summary(kappa.array2)

colnames(minegn.obs.array)<-treesize.array
colnames(minegn.obs1.array)<-treesize.array
colnames(minegn.obs2.array)<-treesize.array


#save.images("fig3.RData")
#load("fig3.rda")
library(gridExtra)
library(ggplot2)
mean.kappa<-apply(log10(kappa.array),2,mean)
mean.kappa1<-apply(log10(kappa.array1),2,mean)
mean.kappa2<-apply(log10(kappa.array2),2,mean)

df <- data.frame(taxa=treesize.array,
                 kappa=c(mean.kappa,mean.kappa1,mean.kappa2), 
                 type=c(rep("birth-death",length(mean.kappa)),rep("coalescent",length(mean.kappa1)),rep("pure birth",length(mean.kappa2))))

legend_title<-"Type"

p1<- ggplot(data=df,aes(x=taxa,y=kappa,group=type,color=type))+
  labs(title="                   kappa vs. taxa (Simulated Data)", y = expression(paste(log[10], " ", kappa)), x="number of taxa")  +
  geom_line(size=1)+
  geom_point(size=2) +
  theme_classic()+
  scale_y_log10(limits=c(3, 10)) 

p1

#minegnobs<-apply(-log10(minegn.obs.array), 2, min)
#minegnobs1<-apply(-log10(minegn.obs1.array), 2, min)
#minegnobs2<-apply(-log10(minegn.obs2.array), 2, min)


#df1 <- data.frame(taxa=treesize.array,
#                 minegn=c(minegnobs,minegnobs1,minegnobs2), 
#                 type=c(rep("birth-death",length(minegnobs)),rep("coalescent",length(minegnobs1)),rep("pure birth",length(minegnobs2))))

#legend_title<-"Type"

#p2<- ggplot(data=df1,aes(x=taxa,y=minegn,group=type,color=type))+
#  labs(title= "Minimum Eigenvalues vs. taxa (Simulated Data))", y = expression(paste(-log[10], " ", min, " ",lambda)), x="number of taxa")  +
#  geom_line(size=1)+
#  geom_point(size=2) +
#  scale_y_log10(limits=c(1, 8))  +
#  theme_classic() 

#p2


#grid.arrange(p1, p2, nrow = 1)

save.image("fig3.rda")