rm(list=ls())
setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationSummary")

rm(list=ls())

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

mu<-0
sigmasq<-1
treesize.array<-c(100,500,800)
simtrait<-50

solve_iid_RMSD<-NULL
for(treesizeIndex in 1:length(treesize.array)){
  treesize<-treesize.array[treesizeIndex]
  for(simtreeIndex in 1: 100){
    print(simtreeIndex)
    brith_lambda=runif(1,0.01,0.1)
    death_mu=runif(1,0,brith_lambda)
    tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
    # shrink.tol<- runif(1, min=1e-15, max=5e-13)
    # C <- tinytipvcv(C=vcv(tree),shrink=shrink.tol) # Bad Tree we use a very bad tree to simulate data, as expected parameter estimate will be bad. so adjust tree will lead a better estimate. need to determine the bound for the tree.
    Yset<-array(0,c(treesize,simtrait))
    C<- vcv(tree)
    one<-array(1,c(dim(C)[1],1))
    for (simtraitIndex in 1:simtrait) {
      Yset[,simtraitIndex]<- mvrnorm(n = 1, mu*one, Sigma=sigmasq*C, tol = 1e-10, empirical = FALSE, EISPACK = FALSE)
      }
    C<-diag(dim(C)[1])
    try(solve_iid_RMSD<-rbind(solve_iid_RMSD,RMSD_sametrait(Yset=Yset, adjC=C, droptipname=NULL, method="solve")))
    }  
}

save.image("rmsdShrinkDiscussion.rda")

boxplot(solve_iid_RMSD)
head(solve_iid_RMSD)
dim(solve_iid_RMSD)

df.iid.plot<-data.frame(RMSD =c(solve_iid_RMSD[1:100,1], solve_iid_RMSD[1:100,2],
                                solve_iid_RMSD[101:200,1], solve_iid_RMSD[101:200,2],
                                solve_iid_RMSD[201:300,1], solve_iid_RMSD[201:300,2]), 
                        Param = c(rep("mu",100),
                                  rep("sigmasq",100),
                                  rep("mu",100),
                                  rep("sigmasq",100),
                                  rep("mu",100),
                                  rep("sigmasq",100)
                                  ),
                        Taxa = c(rep(100,200),
                                 rep(500,200),
                                 rep(800,200)
                                 )
                        )


df.iid.plot$Taxa <- factor(df.iid.plot$Taxa, levels = c("100","500","800"))
p.iid<-ggplot(df.iid.plot,aes(x=Taxa,y=RMSD)) + geom_boxplot(aes(fill=Param)) + scale_y_continuous(trans = "log10")
p.iid<- p.iid + ggtitle(label = paste("Shrinkage Method:  delta = 1"))#, subtitle = paste("θ = ", mu , ", n = ", size ,sep=""))
p.iid<- p.iid + theme(
  plot.title = element_text(color = "black", size = 20, face = "bold"))+ theme(plot.title = element_text(hjust = 0.5))
p.iid
