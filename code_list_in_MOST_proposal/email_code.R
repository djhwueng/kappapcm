#Brian real datasets
rm(list=ls())
library(ape)
myenv<-new.env()
load("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/opentree_chronograms.rda",env=myenv)
ls(envir=myenv)
myenv$datelife.cache$trees[1]$`Arbabi, Tayebeh, Javier Gonzalez, Michael Wink. 2014. A re-evaluation of phylogenetic relationships within reed warblers (Aves: Acrocephalidae) based on eight molecular loci and ISSR profiles. Molecular Phylogenetics and Evolution 78: 304-313.`$tip.label

taxa.kappa<-array(0,c(126,2))

for(index in 1:126){
  or
  tree<-myenv$datelife.cache$trees[[index]]
  taxa.kappa[index,1]<-length(tree$tip.label)
  taxa.kappa[index,2]<-kappa(vcv(tree))
  }

load("/Users/djhwueng/Downloads/opentree_chronograms.rda")
ls(datelife.cache)
datelife.cache$authors[[93]]
datelife.cache$studies[[73]]
datelife.cache$curators[[73]]

datelife.cache$trees[[73]]$edge.length

size.array<-c(0,c(126))
for(i in 1:126){
size.array[i]<-(length(datelife.cache$trees[[i]]$tip.label))
}

summary(size.array)

#Dean shrinkage estimator
rm(list=ls())
library(TreeSim)
library(Matrix)
library(phytools)
#Theiler's Shrinkage Estimator (2012)
TH.cov.shrink<-function(S,n){
  T<-diag(diag(S))
  alpha<-10^(optimize(avenegLOOL,c(-8,0),S,T,n)$minimum)
  sig.pk<-(1-alpha)*S + alpha*T
  sig.pk
}

TH.shrink.alpha<-function(S,n){
  T<-diag(diag(S))
  alpha<-10^(optimize(avenegLOOL,c(-8,0),S,T,n)$minimum)
  return(alpha)
}

#Theiler 2012 optimization to find alpha
avenegLOOL<-function(logalpha,S,T,n){ 
  alpha<-10^logalpha ;  beta<-(1-alpha)/(n-1)
  G<-beta*n*S + alpha*T
  elu<-expand(lu(G))
  L<- elu$P%*%elu$L ; U<-elu$U
  logGdet<-sum(log(abs(diag(U))))
  ro<-sum(diag(solve(G)%*%S)) 
  if (beta*ro>1){
    val<-Inf
  }else{
    val<-log(1-beta*ro) + ro/(1-beta*ro) + logGdet
  }
  return(val)
}


trees <- sim.bd.taxa.age(100, 100, lambda=0.4, mu=0.1, frac=.5, age=log(100)/0.3)
plot(trees[[1]])
vcv.obs<-lapply(1:100, function(x) vcv(trees[[x]]))
kappa.obs<-unlist(lapply(1:100, function(x) kappa(vcv.obs[[x]],exact=TRUE)))
head(kappa.obs)
#Theiler 2012 shrinkage method
vcv.shrunk<-lapply(1:100, function(x) TH.cov.shrink(vcv.obs[[x]],100))
dim(vcv.shrunk[[1]])
kappa.shrunk<-unlist(lapply(1:100, function(x) kappa(vcv.shrunk[[x]],exact=TRUE)))

est.alpha<-unlist(lapply(1:100, function(x) TH.shrink.alpha(vcv.obs[[x]],100)))

plot(est.alpha)#the mle for shrinkage estimate alpha
print(summary(est.alpha))

plot(log(kappa.obs),log(kappa.shrunk),xlab=expression(paste("log ",kappa,(C),sep=""))
     , ylab=expression(paste("log ",kappa,(C[alpha]),sep="")), main="Shrunk Trees vs. Raw Trees"
     ) #Really works well.
axis(1,col="blue",col.axis="blue")
axis(2,col="red",col.axis="red")

head(log(kappa.obs))

head(log(kappa.shrunk))

#pulls branches off 1:1 line
#I want to visualize the before and after as trees
plot(vcv.obs[[1]],vcv.shrunk[[1]])

abline(0,1)

library(phangorn)
#par(mfrow=c(1,2))
for(index in 1:100){
  #index<-1
  print(index)
  shrunk.tree<-upgma(2*(max(vcv.shrunk[[index]])-vcv.shrunk[[index]])) #shall be doubled as it is the distance
  plot(cophylo(trees[[index]], shrunk.tree))
  Sys.sleep(5)
  }

#Brian drop tip
rm(list=ls())
library(ape)
library(TreeSim)

DropShortest <- function(phy) {
  tip.lengths <- phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  phy <- drop.tip(phy, phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]])
  return(phy)
}

getShortest<-function(phy){
  tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  return(which.min(tip.lengths))
  }

DropRandom <- function(phy) {
  phy <- drop.tip(phy, sample.int(Ntip(phy),1))
  return(phy)		
}

CompareKappas <- function(phy) {
  results <- c(kappa(vcv(phy), exact=TRUE), kappa(vcv(DropShortest(phy)), exact=TRUE), kappa(vcv(DropRandom(phy)), exact=TRUE))
  names(results) <- c("Original", "DropShortest", "DropRandom")
  return(results)
}

GetMinEigen<-function(phy){
  evalues<-eigen(vcv(phy))$values
  return(which(evalues==min(evalues))==getShortest(phy))
  }

trees <- sim.bd.taxa.age(100, 100, lambda=0.4, mu=0.1, frac=.5, age=log(100)/0.3)

MinEvalueDropMinTip<-sapply(trees,GetMinEigen)




comparisons <- sapply(trees, CompareKappas)

plot(x=range(comparisons), y=range(comparisons), xlab="original kappa", ylab="kappa: drop 1 tip", type="n", bty="n", log="xy", main="Drop shortest tip vs. Original")
for (i in sequence(dim(comparisons)[2])) {
  lines(rep(comparisons[1,i],2), 	comparisons[2:3,i])
  points(comparisons[1,i], comparisons[2,i], pch=16, col="red")
  points(comparisons[1,i], comparisons[3,i], pch=22, col="blue")
  
}





TestUltrametric <- function(x) {
  ultrametric=TRUE
  for (i in sequence(dim(x)[1])) {
    for (j in sequence(dim(x)[2])) {
      if(x[i,j] < min(c(x[i,], x[,j]))) {
        ultrametric=FALSE
      }
      if (i!=j) {
        if(x[i,i] < x[i,j]) {
          ultrametric=FALSE
        }
      }
    }
  }
  return(ultrametric)
}


for(index in 1:100){
print(TestUltrametric(vcv(trees[[index]])))
}


vcv(trees[[index]])[1:10,1:10]


#Mike add noise to singular value/eigenvalue
rm(list=ls())
library(ape)
library(phytools)
library(MASS)

phy <- pbtree(n=6)
?pbtree
plot(phy)
C <- vcv.phylo(phy)
C
svdC <- svd(C)
U <- svdC$u
D <- diag(svdC$d)
N <- nrow(U)

U%*%D%*%t(U) # same as C

noise <- mvrnorm(N, mu = rep(0,N), Sigma = diag(0.01,N)) # 
noise
U
#small error
UU <- U + noise
UU[which(U==0)] = 0
UU[which(U==1)] = 1

newC <- UU%*%D%*%t(UU)

newC
C

newd <- svd(newC)$d
d <- svdC$d

# distribution of eignevalues should be more palatable, but might have to redo
# also, might not be better if a fix is not needed
newd
d





#Brian  to get distance, we shall double the value.
library(phytools)
library(phangorn)

phy <- pbtree(n=30)
C <- vcv.phylo(phy)
newC <- C
newD <- 2* (max(newC)-newC)

shrunk.tree<-upgma(newD)
#plot(phy)
#plot(shrunk.tree)
plot(cophylo(phy, shrunk.tree))
sort(branching.times(shrunk.tree)) - sort(branching.times(phy))



#shortest tip lengths  = smallest eigen values
rm(list=ls())
library(phytools)
library(phangorn)
phy <- pbtree(n=100)
tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
C <- vcv.phylo(phy)
print(min(tip.lengths)-min(eigen(C)$values))


#omega: condition of the tree
rm(list=ls())
library(TreeSim)
library(ape)

#quadratic<-function(phy){
#  G<-vcv(phy)
#  t(y)%*%solve(G,y)
#  }

kappaget<-function(phy){
  kappa(vcv(phy))
  }

DropShortest <- function(phy) {
  tip.lengths <- phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  phy <- drop.tip(phy, phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]])
  return(phy)
}

sims<-1000
#y<-rnorm(100)
trees <- sim.bd.taxa.age(100, sims, lambda=0.4, mu=0.1, frac=.5, age=log(100)/0.3)
#plot(sapply(trees,kappaget),sapply(trees,quadratic))
bad.trees<-subset(trees, sapply(trees,kappaget)>1e04)
#plot(sapply(illtrees,kappaget),sapply(illtrees,quadratic))
good.tree<-subset(trees, sapply(trees,kappaget)<500)
a.bad.tree<-subset(trees, sapply(trees,kappaget)>1500000)

#t(y)%*%solve(vcv(agoodtree[[1]]),y)
#t(y)%*%solve(vcv(abadtree[[1]]),y)
#min(eigen(vcv(abadtree[[1]]))$values)
#mean(agoodtree[[1]]$edge.length)/min(agoodtree[[1]]$edge.length[which(agoodtree[[1]]$edge[,2]<=Ntip(agoodtree[[1]]))])



drop.a.tip.bad.tree<-DropShortest(a.bad.tree[[1]])

omega.goodtree<- min(eigen(vcv(good.tree[[1]]))$values) /diag(vcv(good.tree[[1]]))[1]
omega.badtree<- min(eigen(vcv(a.bad.tree[[1]]))$values) /diag(vcv(a.bad.tree[[1]]))[1]
omega.drop.badtree<-min(eigen(vcv(drop.a.tip.bad.tree))$values) /diag(vcv(drop.a.tip.bad.tree))[1]

omega.goodtree
omega.badtree
omega.drop.badtree



par(mfrow=c(1,3))
plot(a.bad.tree[[1]],main=paste("bad tree, kappa = ",
                                round(kappa(vcv(a.bad.tree[[1]]))),sep=""))
plot(good.tree[[1]],main=paste("good tree, kappa = ",
                               round(kappa(vcv(good.tree[[1]]))),sep=""))

plot(drop.a.tip.bad.tree,main=paste("drop, kappa = ",
                               round(kappa(vcv(drop.a.tip.bad.tree))),sep=""))


#condition number upper bound test
test.kappa<-function(eps){
  M<-matrix(c(1,-1,-1,1+eps),ncol=2)
  y<-matrix(c(2,3),ncol=1)
  return( list(cond=kappa(M),quadratic=c(array(t(y)%*%solve(M,y)))))
  }

test.kappa(0.000001)



test<-rpois(5,10)
test<-sort(test)
print(test[5]/test[1]>test[4]/test[2])



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


sims<-10 #both replicates of  tree and trait
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

save.images("comparison.RData")


rm(list=ls())
load("comparison.RData")
kappaget<-function(phy){
  kappa(vcv(phy))
}
plot( log(sapply(trees,kappaget)),apply(results$raw,2,rmse,true=sigsq))
plot( log(sapply(trees,kappaget)),apply(results$drop,2,rmse,true=sigsq))
plot( log(sapply(trees,kappaget)),apply(results$shrink,2,rmse,true=sigsq))

summary(lm(apply(results$raw,2,rmse,true=sigsq) ~ log(sapply(trees,kappaget))))

summary(lm(apply(results$drop,2,rmse,true=sigsq) ~ log(sapply(trees,kappaget))))


for(i in 1:50){
   print(fitC.shrink[[i]]$opt$lnL)
}
treeIndex

kappa(vcv(trees[[treeIndex]]))

phytree<-sim.bd.taxa.age(n=100, numbsim=1, lambda=1, mu=0.5, frac = 0.5, age=50, mrca = TRUE)[[1]]
plot(phytree)
omega(phytree)
kappa(vcv(phytree))


fitContinuous(phy=network$phy, dat=means.modified, model="BM")

plot(network$phy)


#kappa increase with N
rm(list=ls())
setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/")
library(TreeSim)
library(ape)
library(phytools)
kappaget<-function(phy){
  kappa(vcv(phy))
  }
treesize.array<-seq(50,850,by=50)
numsim<-100
kappa.array<-array(0,c(numsim,length(treesize.array)))#sim.bd.taxa.age method
kappa.array1<-array(0,c(numsim,length(treesize.array))) #rtree method
kappa.array2<-array(0,c(numsim,length(treesize.array))) #pbtree method
#?sim.bd.taxa.age
#?pbtree
#pbtree(b=1,d=0,n=10)
for(treesizeIndex in 1:length(treesize.array)){
  #tree
  trees <- sim.bd.taxa.age(n=treesize.array[treesizeIndex], numbsim=numsim, lambda=0.4, mu=0.1, frac=.5, age=log(100)/0.3)
  vcv.obs<-lapply(1:numsim, function(x) vcv(trees[[x]]))
  kappa.obs<-unlist(lapply(1:numsim, function(x) kappa(vcv.obs[[x]],exact=TRUE))) 
  kappa.array[,treesizeIndex]<-kappa.obs
  #tree1
  trees1<-    rep(rtree(treesize.array[treesizeIndex]),numsim)
  vcv.obs1<-lapply(1:numsim, function(x) vcv(trees1[[x]]))
  kappa.obs1<-unlist(lapply(1:numsim, function(x) kappa(vcv.obs1[[x]],exact=TRUE))) 
  kappa.array1[,treesizeIndex]<-kappa.obs1
  #tree 2
  trees2<-    rep(pbtree(n=treesize.array[treesizeIndex]),numsim)
  vcv.obs2<-lapply(1:numsim, function(x) vcv(trees2[[x]]))
  kappa.obs2<-unlist(lapply(1:numsim, function(x) kappa(vcv.obs2[[x]],exact=TRUE))) 
  kappa.array2[,treesizeIndex]<-kappa.obs2
  }
#rm(list=ls())
log(kappa.array)
kappa_mean1 <-apply(kappa.array,2,mean)
head(log(kappa.array))
log(kappa.array1)
kappa_mean2<-apply(kappa.array1,2,mean)
log(kappa.array2)
kappa_mean3<-apply(kappa.array2,2,mean)
plot(kappa_mean3)
log.kappa.array<-
data.frame(log(kappa.array))
colnames(log.kappa.array)<-treesize.array
#require(ggplot2)
#require(reshape2)
#log.kappa.array<-melt(log.kappa.array)
#log.kappa.array
#p<-ggplot(data=log.kappa.array,aes(x=variable,y=value)) + geom_boxplot(aes(fill=variable))
#p+labs(title="kappa vs. taxa")+labs(x="taxa")+labs(y="log kapa")

med.kappa<-(apply(log(kappa.array),2,median))
med.kappa1<-(apply(log(kappa.array1),2,median))
med.kappa2<-(apply(log(kappa.array2),2,median))

plot(apply(log(kappa.array),2,median),type="b", ylim=c(min(med.kappa,med.kappa2,med.kappa1), max(med.kappa,med.kappa2,med.kappa1)), xlab="number of taxa", ylab= "log kappa", xaxt='n', main="kappa vs. taxa" )

#plot(apply(log(kappa.array),2,median),type="b", ylim=c(min(med.kappa,med.kappa1), max(med.kappa,med.kappa1)), xlab="number of taxa", ylab= "log kappa", xaxt='n', main="kappa vs. taxa", col="red" )
axis(1,at=1:length(treesize.arra),treesize.array)
#lines(1:length(treesize.array),med.kappa1,type="b",pch=2, col="blue")
#lines(1:length(treesize.array),med.kappa2,type="b",pch=3)
#pbtree is not so stable for the 
#?legend
#legend(1,10,"sim.bd",pch=1,title="method")
#legend(, yrange[2], 1:ntrees, cex=0.8, col=colors,
 #      pch=plotchar, lty=linetype, title="Tree")
mean.kappa<-(apply(log(kappa.array),2,mean))
mean.kappa1<-(apply(log(kappa.array1),2,mean))
mean.kappa2<-(apply(log(kappa.array2),2,mean))


save.image("kappa_taxa_sim.RData")

setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/")
load("kappa_taxa_sim.RData")

plot(apply(log(kappa.array),2,median),type="b", ylim=c(min(mean.kappa,mean.kappa2,mean.kappa1), max(mean.kappa,mean.kappa2,mean.kappa1)), xlab="number of taxa", ylab= "log kappa", xaxt='n', main="kappa vs. taxa (Simulated Data)" )

#plot(apply(log(kappa.array),2,median),type="b", ylim=c(min(med.kappa,med.kappa1), max(med.kappa,med.kappa1)), xlab="number of taxa", ylab= "log kappa", xaxt='n', main="kappa vs. taxa", col="red" )
axis(1,at=1:length(treesize.array),treesize.array)
lines(1:length(treesize.array),mean.kappa1,type="b",pch=2, col="blue")
lines(1:length(treesize.array),mean.kappa2,type 


mean.kappa_real.data<-c(6.452163, 7.180684,  9.225768,  9.575948,  8.350919,  9.746075, 10.968950, 12.639862, NA, 14.108984, 11.066064, NA, 12.388033, NA, 13.332289, NA, 11.619834)



#Dean's email that the PB is unstable than the real tree.
plot(mean.kappa2,type="b", ylim=c(min(med.kappa,med.kappa2), max(med.kappa,med.kappa2)), xlab="number of taxa", ylab= "log kappa", xaxt='n', main="kappa vs. taxa", col="red" )
lines(1:length(treesize.array),mean.kappa_real.data,type="b",pch=3, col="purple")

axis(1,at=1:length(treesize.array),treesize.array)



#Brian: real datasets kappa
rm(list=ls())
library(ape)
setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/")
myenv<-new.env()
load("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/opentree_chronograms.rda",env=myenv)
ls(envir=myenv)
#myenv$datelife.cache$trees[1]$`Arbabi, Tayebeh, Javier Gonzalez, Michael Wink. 2014. A re-evaluation of phylogenetic relationships within reed warblers (Aves: Acrocephalidae) based on eight molecular loci and ISSR profiles. Molecular Phylogenetics and Evolution 78: 304-313.`$tip.label
kappa.array<-array(0,c(130))
size.array<-array(0,c(130))

for(treeIndex in 1:130){
  print(treeIndex)
#  if(treeIndex !=34){
  vcv.matrix<-vcv(myenv$datelife.cache$trees[[treeIndex]])
  kappa.array[treeIndex]  <- kappa(vcv.matrix, exact=TRUE)
  size.array[treeIndex] <- length(myenv$datelife.cache$trees[[treeIndex]]$tip.label)
#  }
  }

kappa.array<-kappa.array[-34]
size.array<-size.array[-34]
save.image("real_data_kappa.RData")
setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/")
load("real_data_kappa.RData")
kappa.array<-kappa.array[-(23:25)]
size.array<-size.array[-(23:25)]
plot(size.array,log(kappa.array),ylab="log kappa", xlab="taxa", main="condition number vs. taxa",lwd=1,lty=1,pch=5)

rm(list=ls())
setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/")
load("real_data_kappa.RData")
in.order<-order(size.array)
plot(size.array[in.order[1:126]],log(kappa.array[in.order][1:126]),xlab="number of taxa",ylab="log kappa",main="kappa vs. taxa (Real Data)",col=4, lwd=2)
in.order<-order(size.array)
kappa.array<-log(matrix(kappa.array[in.order], ncol=1))
size.array<-matrix(size.array[in.order], ncol=1)
output<-cbind(size.array,kappa.array)
colnames(output)<-c("size","log_kappa")
output<-data.frame(output)
output<-output[1:(dim(output)[1]-3),]
output$bin<-cut(output$size,17, labels=F)
agg.output<-aggregate(log_kappa~bin,data=output ,mean)
plot(agg.output$bin,agg.output$log_kappa, type="b", xlab="taxa",ylab="log kappa",main="Real Data: kappa vs. taxa", xaxt='n')
length(output$bin)

agg.output$log_kapp

agg.output<-c(6.452163, 7.180684,  9.225768,  9.575948,  8.350919,  9.746075, 10.968950, 12.639862, NA, 14.108984, 11.066064, NA, 12.388033, NA, 13.332289, NA, 11.619834)

length(agg.output)

axis(1,at=1:17,c(50,100,150,200,250,300,350,400,NA,450,500,NA,550,NA,650,NA,800))



illcond<-function(eps,n=n,Y=Y){
  mtx<-diag(1,n+1)
  C<-matrix(1,ncol=n,nrow=n)
  diag(C)<-diag(C)+eps
  mtx[1:dim(C)[1],1:dim(C)[1]]<-C
  print(solve(mtx,Y))
  return(kappa(mtx))
}
n<-2
Y<-matrix(rnorm(n+1),ncol=1)
eps<-0.000001
illcond(eps,n=n,Y=Y)

eps<-seq(0.001,0.01,0.001)
sapply(eps,illcond,n=n,Y=Y)


library(ape)
tree<-rcoal(3)
plot(tree)
vcv(tree)

# when the search starts, the vcv changes each time the parameters vary and kappa change accordingly. What to do #before analysis: drop tip, strech length, lengthed tips  
#during analysis: 
#after analysis:




#### Length Stretching

rm(list=ls())
#----------------------------------------------
SegM<-function(M){#find branch length segments
  segment<-NULL
  sortM<-sort(M)
  while(sum(sortM!=0)){ #from tip to bottom 
    m<-max(sortM)        
    sortM[which(sortM==m)]<-0
    sm<-max(sortM)
    segment<-cbind(segment,m-sm)# t_d,...,t_2,t_1
  }
  return(rev(segment)) #t_1,t_2,...,t_d from root to tips.
}

PtbLen<-function(sgmt){ #varying time segments
  k<-ceiling(1/min(sgmt)) 
  applen<-array(0,c(length(sgmt)))
  for (i in 1:length(sgmt) ){
    applen[i]<-rgamma(1,shape=k*sgmt[i],rate=1) 
  }
  return(applen/sum(applen)) #Dirichlet samples
}


PtbLen_beta<-function(sgmt){ #varying time segments under beta
  new_sgmt<-array(0,length(sgmt))
  new_sgmt[1]<-rbeta(n=1,shape1 = sgmt[1],shape2 = sum(sgmt[-1]))
  for(Index in 2:(length(sgmt)-1)){
    s_Index<-rbeta(n=1,shape1 = sgmt[Index],shape2 = sum(sgmt[-(1:Index)]))
    new_sgmt[Index]<-(1-sum(sgmt[1:(Index-1)]))*s_Index
  }
  new_sgmt[length(sgmt)] <- 1 - sum(new_sgmt[-length(new_sgmt)]) 
  return(new_sgmt)
}


AccSum<-function(sgmt){ #accumulate segments
  accsum<-array(0,c(length(sgmt))) 
  tmp<-0
  for(i in 1:length(sgmt)){
    tmp<-tmp+sgmt[i]   
    accsum[i]<-tmp
  }
  return(accsum)
} #g1=t1,g2=t1+t2,...,gd=t1+t2+...+td


VOR<-function(accsum, mtx=mtx){ #obtain the matrix
  accsum<-rev(accsum) #gd, ... ,g2 ,g1
  outmtx<-array(0,c(dim(mtx)))
  M<-mtx
  count<-0
  while(sum(M)!=0){
    count<-count+1
    ptb<-max(M)
    positions<-which(mtx==ptb)
    outmtx[positions]<-accsum[count]   
    M[positions]<-0  
  }
  return(outmtx)
}  

PerturbTree<-function(phy){ #get the new tree
  mtx<-vcv(phy)
  height<-max(mtx)
  mtx<-mtx/height #scale tree to make the Dirichlet distribution work
  sgmt<- SegM(mtx) #get segments
  ptblen<-PtbLen(sgmt) #perturb segments
  accsum<-AccSum(ptblen)
  M<-VOR(accsum,mtx=mtx)
  rownames(M)<-rownames(mtx)
  colnames(M)<-colnames(mtx)
  M<-height*M # put the tree height back for the perturb tree
  phy<-upgma(2*(max(M)-M))
  return(phy)
}

kappaget<-function(phy){
  kappa(vcv(phy))
}

#library(phyclust)
#library(geiger)
library(TreeSim)
library(ape)
#library(Matrix)
#library(phytools)
library(phangorn)
#library(geiger)
#library(MASS)


treesize<-20
repeat_n<-500
numbsim<-1 #simulate 1 tree, treat it as the true tree.
trees<-sim.bd.taxa.age(n=treesize,numbsim=numbsim,lambda=1,mu=0.5, frac=0.5, age= (log(treesize)/(1-0.5)))
ptrees<-PerturbTree(trees[[1]])
pKappaset<-replicate(n=repeat_n,kappaget(PerturbTree(trees[[1]]))) # perturb the tree by stretching the lengths. obtain 100 kappa.  
sum(log(kappa(vcv(trees[[1]])))<log(pKappaset))/repeat_n #imporve kappa


pKappaset<-replicate(n=repeat_n,kappaget(PerturbTree(trees[[1]]))) # perturb the tree by stretching the lengths. obtain 100 kappa.  
pKappaset_beta<-replicate(n=repeat_n,kappaget(PerturbTree_beta(trees[[1]]))) # perturb the tree by stretching the lengths. obtain 100 kappa.  

par(mfrow=c(1,2))
plot(log(pKappaset),main="Length Stretch: Gamma",ylab="log kappa", xlab="replicates index")
abline(h=log(kappa(vcv(trees[[1]]))),col="red")
plot(log(pKappaset_beta),main="Length Stretcg: Beta",ylab="log kappa", xlab="replicates index")
abline(h=log(kappa(vcv(trees[[1]]))),col="red")

par(mfrow=c(1,3))
plot(trees[[1]],main="raw tree")
plot(PerturbTree(trees[[1]]),main="new tree:gamma")
plot(PerturbTree_beta(trees[[1]]),main="new tree:beta")


