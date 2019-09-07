setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationRun3/")
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

############## shrink method
avenegLOOL<-function(logdelta,C,T,n){
  delta<-10^logdelta
  beta <-(1-delta)/(n-1)
  S_delta<-n*beta*C+delta*T #formula (3) in
#  elu<-expand(lu(S_delta))
  elu<-lu(S_delta)
  #elu<-lu(S_delta)
  #L<-elu$P%*%elu$L
  L<-elu$L
  U<-elu$U
  logSdet<-sum(log(abs(diag(U))))
  ro<-sum(diag(solve(S_delta)%*%C))
  if(beta*ro>1){
    val<-Inf
  }else{
    val<-log(1-beta*ro)+ro/(1-beta*ro)+logSdet
  }
  return(val)
}


TH.cov.shrink<-function(S,n){
  T<-diag(diag(S))
  delta<-10^(optimize(avenegLOOL,c(-10,0),S,T,n)$minimum)
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


Shrink.C<-function(C){
  phy.vcv<-C
  vcv.shrunk<-TH.cov.shrink(phy.vcv,dim(phy.vcv)[1])
  return(vcv.shrunk/max(vcv.shrunk))
}


############# Drop Tip Method


DropShortest<-function(phy,traitset){
  tip.lengths <-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  drop.tip.index<-phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lenths)[1]]
  drop.tip.names<-names(tip.lengths[drop.tip.index])
  phy<-drop.tip(phy,drop.tip.index)
  traitset<-traitset[-which(rownames(traitset)==drop.tip.names),]
  return(list(phy=phy,traitset=traitset))
}


DropShortesttip<-function(phy){
  tip.lengths <-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  drop.tip.index<-phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]]
  drop.tip.names<-names(tip.lengths[drop.tip.index])
  phy<-drop.tip(phy,drop.tip.index)
  return(list(phy=phy,drop.tip.names=drop.tip.names))
}


##############Stretch Length Method

SegM<-function(M){#find branch length segments
  segment<-NULL
  sortM<-sort(M)
  while(sum(sortM!=0)){#from tip to bottom
    m<-max(sortM)
    sortM[which(sortM==m)]<-0
    sm<-max(sortM)
    segment<-cbind(segment,m-sm)#t_d,...,t_2,t_1
  }
  return(rev(segment))#t_1,t_2,...,t_d from root to tips
}

PtbLen_gamma<-function(sgmt){#varying time segments
  k<-ceiling(1/min(sgmt))
  applen<-array(0,c(length(sgmt)))
  for (i in 1:length(sgmt)){
    applen[i]<-rgamma(1,shape=k*sgmt[i],rate=1)
  }
  return(applen/sum(applen)) #Dirichlet samples
}

AccSum<-function(sgmt){#accumulate segments
  accsum<-array(0,c(length(sgmt)))
  tmp<-0
  for(i in 1:length(sgmt)){
    tmp<-tmp+sgmt[i]
    accsum[i]<-tmp
  }
  return(accsum)
} #g1=t1, g2=t1+t2,...,gd=t1+t2+...+td

VOR<-function(accsum,mtx=mtx){ #obtain the matrix
  accsum<-rev(accsum)#gd,...,g2,g1
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

PerturbTree<-function(phy,method=method){#get the new tree
  mtx<-vcv(phy)
  height<-max(mtx)
  mtx<-mtx/height #scale tree to make the Dirichlet distribution work
  sgmt<-SegM(mtx) #get segmetns
  if(method=="gamma"){
    ptblen<-PtbLen_gamma(sgmt) #perturb segments
  }
  if(method=="beta"){
    ptblen<-PtbLen_beta(sgmt) #perturb segments
  }
  accsum <-AccSum(ptblen)
  M<-VOR(accsum,mtx=mtx)
  rownames(M)<-rownames(mtx)
  colnames(M)<-colnames(mtx)
  M<-height*M #put the tree height back for the perturb tree
  phy<-upgma(2*(max(M)-M))
  return(phy)
}


PerturbC<-function(C){#get the new tree
  mtx<-C
  height<-max(mtx)
  mtx<-mtx/height #scale tree to make the Dirichlet distribution work
  sgmt<-SegM(mtx) #get segmetns
  ptblen<-PtbLen_gamma(sgmt) #perturb segments
  #print(PtbLen_gamma(sgmt))
  accsum <-AccSum(ptblen)
  M<-VOR(accsum,mtx=mtx)
  rownames(M)<-rownames(mtx)
  colnames(M)<-colnames(mtx)
  M<-height*M #put the tree height back for the perturb tree
  return(M/max(M))
}


DropShortesttipName<-function(phy){
  tip.lengths <-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  names(tip.lengths)<-phy$tip.label
  drop.tip.index<-phy$edge[which(phy$edge[,2]<=Ntip(phy)),2][which.min(tip.lengths)[1]]
  drop.tip.names<-names(tip.lengths[drop.tip.index])
  return(drop.tip.names)
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


###########root-mean-square deviation
RMSD<-function(true.model.param=true.model.param,sims=sims,C=C,adjC=adjC,droptipname=NULL,method=method){
  true.model.param[1]<-mu
  true.model.param[2]<-sigmasq
  param.array<-array(0,c(sims,2))
  one<-array(1,c(dim(C)[1],1))
  for (i in 1:sims) {
    y<- mvrnorm(n = 1, mu*one, Sigma=sigmasq*C, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    names(y)<-rownames(C)
    if(!is.null(droptipname)){
      y<-y[-(names(y)==droptipname)]
    }
    param.array[i,]<- MLEfcn(adjC=adjC,Y=y,method=method)
  }
  colnames(param.array)<-c("theta","sigmasq")

  rmsd_mu<-  sqrt(mean((param.array[,1]-mu) ,na.rm=TRUE)^2)    #sqrt(mean((param.array[,1]-mu)^2))/mean(param.array[,1],na.rm=TRUE)
  rmsd_sigmasq <- sqrt(mean(param.array[,2]-sigmasq, na.rm=TRUE)^2)#sqrt(mean((param.array[,2]-sigmasq)^2))/mean(param.array[,2],na.rm=TRUE)
  return( c(rmsd_mu,rmsd_sigmasq))
}

svd_decom<-function(M){#returns M^{-1}.
  decom<-svd(M)
  val<-decom$d
  solveD<-diag(1/val)
  return(decom$v%*%solveD%*%t(decom$u))
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


RMSD_geiger<-function(Yset=Yset,tree=tree){
  sims<-dim(Yset)[2]
  param.array<-array(0,c(sims,2))
  for (i in 1:sims) {
    y<-Yset[,i]
    names(y)<-tree$tip.label
    gfit<-fitContinuous(phy=treeupgma,dat=y, model=c("BM"))
    param.array[i,]<- c(gfit$opt$z0,gfit$opt$sigsq)
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
solCset<-NULL
unsolCset<-NULL
soltreeset<-NULL
unsoltreeset<-NULL

countsol<-0
countunsol<-0

treesize<-100 # 100, 500, 800
simtrait<-50 # 50 traits
Yset<-array(0,c(treesize,simtrait))
mu<-0 # or mu <- 10
sigmasq<-1# sigmasq<- 5

savename<-paste("size",treesize,"mu",mu,"sigmasq",sigmasq,".rda",sep="")
savename

for(treeIndex in 1:simsTree){
    print(treeIndex)
    solveC<-NULL
    brith_lambda=runif(1,0.01,0.1)
    death_mu=runif(1,0,brith_lambda)
    tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
    shrink.tol<- runif(1, min=1e-15, max=5e-13)
    C <-  tinytipvcv(C=vcv(tree),shrink=shrink.tol) # Bad Tree we use a very bad tree to simulate data, as expected parameter estimate will be bad. so adjust tree will lead a better estimate. need to determine the bound for the tree.
    kappaC<-kappa(C)

    try(solveC<-solve(C))
    if (is.null(solveC)){
      countunsol <- countunsol + 1
      unsoltreeset[[countunsol]]<-tree
      unsolCset[[countunsol]] <- C/max(C)
      }else{
        countsol <- countsol + 1
        soltreeset[[countsol]]<-tree
        solCset[[countsol]]<-C/max(C)
    }
  }

print(countsol)
print(countunsol)



# Do RMSD for solvable tree
solve_raw_RMSD<-NULL
solve_shrink_RMSD<-NULL
solve_stretch_RMSD<-NULL
solve_droptip_RMSD<-NULL
solve_prun_RMSD<-NULL


for(solIndex in 1: countsol){
#  solIndex<-2
  print(solIndex)
  C<-solCset[[solIndex]]
  one<-array(1,c(dim(C)[1],1))
  for (simtraitIndex in 1:simtrait) {
    Yset[,simtraitIndex]<- mvrnorm(n = 1, mu*one, Sigma=sigmasq*C, tol = 1e-10, empirical = FALSE, EISPACK = FALSE)
  }
  #1 solvable tree
 # try(print(RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  try(solve_raw_RMSD<-rbind(solve_raw_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  #2. Shrinkage
  C<-Shrink.C(solCset[[solIndex]])
#  print(RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve"))
  try(solve_shrink_RMSD<-rbind(solve_shrink_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  #3 Stretching
  solveC<-NULL
  attemp<-0
  while(is.null(solveC) & attemp<100){
    C<-solCset[[solIndex]]
    C<-PerturbC(C)
    attemp=attemp+1
  try(solveC<-solve(C))
  }
#  try(print(RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  try(solve_stretch_RMSD<-rbind(solve_stretch_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))

  #5 pruning algoruthm
  C<-solCset[[solIndex]]
  D<- 2*(diag(C)[1]-C)
  treeupgma<-upgma(D)
#  try(print(RMSD_geiger(Yset=Yset,tree=treeupgma)))
  try(solve_prun_RMSD<-rbind(solve_prun_RMSD,RMSD_geiger(Yset=Yset,tree=treeupgma)))

  #4  Droptip
  droptipname<- DropShortesttipName(soltreeset[[solIndex]])
  delc<-which(rownames(solCset[[solIndex]])==droptipname)
  C<- solCset[[solIndex]][-delc,-delc]
  one<-array(1,c(dim(C)[1],1))
  dYset<-Yset[-delc,]
#  try(print(RMSD_sametrait(Yset=dYset, adjC=C, droptipname=NULL, method="solve")))
  try(solve_droptip_RMSD<-rbind(solve_droptip_RMSD,RMSD_sametrait(Yset=dYset, adjC=C, droptipname=NULL, method="solve")))

save.image(savename)
}

# Do RMSD for unsolvable tree
unsolve_shrink_RMSD<-NULL
unsolve_stretch_RMSD<-NULL
unsolve_droptip_RMSD<-NULL
unsolve_prun_RMSD<-NULL

for(solIndex in 1: countunsol){
  #solIndex<-
  print(solIndex)
  C<-unsolCset[[solIndex]]
  one<-array(1,c(dim(C)[1],1))
  for (simtraitIndex in 1:simtrait) {
    Yset[,simtraitIndex]<- mvrnorm(n = 1, mu*one, Sigma=sigmasq*C, tol = 1e-10, empirical = FALSE, EISPACK = FALSE)
  }
  #1 solvable tree
  #try(print(RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  #try(unsolve_raw_RMSD<-rbind(unsolve_raw_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  #2. Shrinkage
  C<-Shrink.C(unsolCset[[solIndex]])
#  print(RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve"))
  try(unsolve_shrink_RMSD<-rbind(unsolve_shrink_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  #3 Stretching
  solveC<-NULL
  attemp<-0
  while(is.null(solveC) & attemp<100){
    C<-unsolCset[[solIndex]]
    C<-PerturbC(C)
    attemp=attemp+1
  try(solveC<-solve(C))
  }
#  try(print(RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
  try(unsolve_stretch_RMSD<-rbind(unsolve_stretch_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))

  #5 pruning algoruthm
  C<-unsolCset[[solIndex]]
  D<- 2*(diag(C)[1]-C)
  treeupgma<-upgma(D)
#  try(print(RMSD_geiger(Yset=Yset,tree=treeupgma)))
  try(unsolve_prun_RMSD<-rbind(unsolve_prun_RMSD,RMSD_geiger(Yset=Yset,tree=treeupgma)))

  #4  Droptip
  rownames(Yset)<-rownames(unsolCset[[solIndex]])
  droptipname<- DropShortesttipName(unsoltreeset[[solIndex]])
  delc<-which(rownames(unsolCset[[solIndex]])==droptipname)
  C<- unsolCset[[solIndex]][-delc,-delc]
  dYset<-Yset[-delc,]
#  try(print(RMSD_sametrait(Yset=dYset, adjC=C, droptipname=NULL, method="solve")))
  try(unsolve_droptip_RMSD<-rbind(unsolve_droptip_RMSD,RMSD_sametrait(Yset=dYset, adjC=C, droptipname=NULL, method="solve")))

save.image(savename)
}
