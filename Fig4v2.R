setwd("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/")
rm(list=ls())
library(TreeSim)
library(Matrix)
library(phytools)

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

PtbLen_beta<-function(sgmt){#varying time segments under beta
  new_sgmt<-array(0,length(sgmt))
  new_sgmt[1]<-rbeta(n=1,shape1=sgmt[1],shape2=sum(sgmt[-1]))
  for(Index in 2:(length(sgmt)-1)){
    s_Index<-rbeta(n=1,shape1=sgmt[Index],shape2=sum(sgmt[-(1:Index)]))
    new_sgmt[Index]<-(1-sum(sgmt[1:(Index-1)]))*s_Index
  }
  new_sgmt[length(sgmt)]<-1-sum(new_sgmt[-length(new_sgmt)])
  return(new_sgmt)
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


# Shirnkage method

############## shrink method
avenegLOOL<-function(logdelta,C,T,n){
  delta<-10^logdelta
  beta <-(1-delta)/(n-1)
  S_delta<-n*beta*C+delta*T #formula (3) in  
  elu<-Matrix::expand(lu(S_delta))
  #elu<-lu(S_delta)
  L<-elu$P%*%elu$L
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


sims<-200
kappa.obs<-array(0,c(sims))
kappa.droptip<-array(0,c(sims))
kappa.gamma<-array(0,c(sims))
kappa.shrunk<-array(0,c(sims))


for(simIndex in 1:sims){
  #simIndex<-1
  print(simIndex)
  brith_lambda=runif(1,0.01,0.1)
  death_mu=runif(1,0,brith_lambda)
  treesize<-round(runif(1,min=300,max=1000))
#  tree<-sim.bd.taxa.age(n= treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,frac=frac,age=log(treesize)/(brith_lambda-death_mu),mrca=TRUE)[[1]]
  tree<-sim.bd.taxa(n=treesize,numbsim=1,lambda=brith_lambda,mu=death_mu,complete = FALSE, stochsampling = TRUE)[[1]]
  
  
#  treesize<-100
#  tree <- sim.bd.taxa.age(100, 1, lambda=0.4, mu=0.1, frac=.5, age=log(100)/0.3)[[1]]
  
  
  kappa.obs[simIndex]<-kappa(vcv(tree))
  kappa.gamma[simIndex]<-kappa(vcv(PerturbTree(phy=tree,method="gamma")))
  kappa.droptip[simIndex]<-kappa(vcv(DropShortesttip(phy=tree)$phy))
  kappa.shrunk[simIndex]<-kappa(vcv(Shrink.tree(phy=tree)))

  print(round(c(treesize,kappa.obs[simIndex],kappa.gamma[simIndex],kappa.droptip[simIndex],kappa.shrunk[simIndex]  ) ))
}

save.image("Giv4v2.RData")


library(ggplot2)
library(gridExtra)

minaxis<-min(log10(kappa.obs),log10(kappa.shrunk))
maxaxis<-max(log10(kappa.obs),log10(kappa.shrunk))
df<-data.frame(logkappaobs=log10(kappa.obs),logkappashrink=log10(kappa.shrunk))
p1<-ggplot(data=df,aes(x=logkappaobs,y=logkappashrink)) + geom_point(size=2) + labs(x=expression(paste("log10 ", kappa, (C), sep="")), y=expression(paste("log ", kappa, (C[1]), sep=""))) + xlim(minaxis,maxaxis) + ylim(minaxis,maxaxis) + theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + geom_hline(yintercept=max(log(kappa.gamma))+0.1,linetype="dashed",color="red")

p1<- p1 + ggtitle("Shrunk Tree")  + geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1)

minaxis<-min(log10(kappa.obs),log10(kappa.gamma))
maxaxis<-max(log10(kappa.obs),log10(kappa.gamma))
df<-data.frame(logkappaobs=log10(kappa.obs),logkappagamma=log10(kappa.gamma))
p2<-ggplot(data=df,aes(x=logkappaobs,y=logkappagamma)) + geom_point(size=2) + labs(x=expression(paste("log10 ", kappa, (C), sep="")), y=expression(paste("log ", kappa, (C[delta]), sep=""))) + xlim(minaxis,maxaxis) + ylim(minaxis,maxaxis) + theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + geom_hline(yintercept=max(log(kappa.gamma))+0.1,linetype="dashed",color="red")

p2 <- p2 + ggtitle("Length Stretching") + geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1)


minaxis<-min(log10(kappa.obs),log10(kappa.droptip))
maxaxis<-max(log10(kappa.obs),log10(kappa.droptip))
df<-data.frame(logkappaobs=log10(kappa.obs),logkappadroptip=log10(kappa.droptip))
p3<-ggplot(data=df,aes(x=logkappaobs,y=logkappadroptip)) + geom_point(size=2) + labs(x=expression(paste("log10  ", kappa, (C), sep="")), y=expression(paste("log ", kappa, (C[drop]), sep=""))) + xlim(minaxis,maxaxis) + ylim(minaxis,maxaxis) + theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + geom_hline(yintercept=max(log(kappa.droptip))+0.1,linetype="dashed",color="red")
p3 <- p3 + ggtitle("Drop Tip") + geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1)

grid.arrange(p1,p2,p3,ncol=3)

