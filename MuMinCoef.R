setwd("~/Dropbox/JournalSubmission/EvolutionaryBioinformatics-kappa/code_MajorRevision/")
rm(list=ls())
library(TreeSim)
library(ape)
library(MuMIn)

GetSummary<-function(phy,lambda=NA,mu=NA){
  log10.kappa.result<-log10(kappa(vcv(phy,exact=TRUE)))
  min.brlen<-min(phy$edge.length)
  max.brlen<-max(phy$edge.length)
  brlen.ratio<-max.brlen/min.brlen
  brlen.var<-var(phy$edge.length)
  brlen.median<-median(phy$edge.length)
  tip.lengths<-phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
  internal.lengths<-phy$edge.length[which(phy$edge[,2]>Ntip(phy))]
  min.tip<-min(tip.lengths)
  max.tip<-max(tip.lengths)
  min.internal<-min(internal.lengths)
  max.internal<-max(internal.lengths)
  return(data.frame( 
    Ntip=Ntip(phy),
    log10.kappa=log10.kappa.result, 
    min.brlen=min.brlen, 
    max.brlen=max.brlen, 
    brlen.ratio=brlen.ratio, 
    brlen.var=brlen.var,
    brlen.median=brlen.median,
    min.tip=min.tip,
    max.tip=max.tip,
    min.internal=min.internal,
    max.internal=max.internal,
    lambda=lambda,
    mu=mu,
    turnover=1/(lambda+mu)
    ))
  }

nsim<-100
nreps=100

all.results <-data.frame()

for(i in sequence(nreps)){
  ntax <-round(runif(1,10,800))
  lambda=runif(1,0.01,0.1)
  mu=runif(1,0,0.1)
  trees<-sim.bd.taxa(n=ntax,numbsim=nsim,lambda=lambda,mu=mu,complete = FALSE, stochsampling = TRUE)
  tree.results<-t(sapply(trees,GetSummary,lambda=lambda,mu=mu))
  if(i==1){
    all.results <- tree.results
  }else{
    all.results <-rbind(all.results,tree.results)
  }
  print(paste("Done",i,"of",nreps))
}

all.results <-data.frame(all.results)
#head(all.results)
for(i in sequence(dim(all.results)[2])){
  all.results[,i]<-unlist(all.results[,i])
}




options(na.action = 'na.fail')
global.model <-lm(log10.kappa ~ ., data=all.results)
dredge.results <- dredge(global.model, m.lim=c(NA,4),trace=1)
dim(dredge.results)
print(dredge.results,abbrev.names=FALSE)
head(dredge.results)


write.csv(dredge.results,file="dredge.results1_10000trees10to800taxa.csv")
system("open dredge.results1_10000trees10to800taxa.csv")
save.image("MuMinCoef.rda")


predictors<-colnames(dredge.results[,2:14])
TF<-!is.na(dredge.results[,predictors])
table(apply(TF,1,sum))
sum(table(apply(TF,1,sum)))

head(dredge.results)


best.model<-lm(log10.kappa ~ Ntip+ lambda  + mu + min.tip , data=all.results)
summ.best<-summary(best.model)
library(xtable)
xtable(t(summ.best$coefficients[,1:2]),digits=3)
