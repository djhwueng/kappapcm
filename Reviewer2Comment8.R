rm(list=ls())
library(ape)
library(ggplot2)
library(gridExtra)
library(TreeSim)

#setwd("~/Dropbox/ChienSanPeng/Thesis/")
setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/Thesis/")
myenv<-new.env()
load("opentree_chronograms.rda",env=myenv)
ls(envir=myenv)
myenv$datelife.cache$dois

taxa.size.array<-c()
for(index in 1:130){
  taxa.size.array<-c(taxa.size.array,length(myenv$datelife.cache$trees[[index]]$tip.label))
}
taxa.size.array
length(myenv$datelife.cache$trees[[34]]$tip.label)

median(taxa.size.array[-c(23,24,25,34)])
#mean(taxa.size.array[-c(23,24,25,34)])

taxa.kappa.array<-array(0,c(130,2))
for(index in 1:130){
  #  print(index)
  if(index != 23 & index != 24 & index != 25 & index != 34){#exclude the tree with 48016 taxa
    tree<-myenv$datelife.cache$trees[[index]]
    taxa.kappa.array[index,1]<-length(tree$tip.label)
    #print(diag(vcv(tree))[1])
    vcvtree<-vcv(tree)
    taxa.kappa.array[index,2]<-kappa(vcvtree/max(vcvtree))
  }
}

taxa.kappa.array<-taxa.kappa.array[-c(23,24,25,34),]

taxa.kappa.array<-data.frame(taxa=taxa.kappa.array[,1],kappa=taxa.kappa.array[,2])

plot(taxa.kappa.array$taxa, log10(taxa.kappa.array$kappa))


library(ggplot2)
ggplot(data=taxa.kappa.array,mapping = aes(x=taxa,y=log10(kappa))) + geom_point( color = "blue") + 
  geom_smooth(method=lm, se=TRUE, color="black") +
  labs(title="          kappa vs. taxa (Real Data)", y = "log kappa", x="number of taxa")  + theme_classic() + geom_smooth(method=lm)

lm(log10(kappa) ~ taxa, data=taxa.kappa.array)

nsim<-1
numdiffage<-10
# Now we do simulate tree case
#1.Compare the same number of tips between real and simulated data.
sim.kappa.array<-taxa.size.array[-c(23,24,25,34)]
sim.taxa.kappa.array<-array(0,c(length(sim.kappa.array),numdiffage,2))
#2.Take into account differences in tree height (i.e., distance from the root to the tips of the tree. 
# For each number of tips, we consider 100 trees in differences tree height
for(sizeIndex in 1:length(sim.kappa.array)){
  print(sizeIndex)
  ntax <-sim.kappa.array[sizeIndex]
  for (ageIndex in 1:numdiffage){
    print(ageIndex)
    lambda=runif(1,0.01,0.1)
    mu=runif(1,0,lambda)
    #frac=runif(1,0.1,1)
    #age=log(ntax)/(lambda-mu)
    # 3. Make sure that all trees are ultrametic (if not, then discuss if there is any difference.)
    # Here we are sure that all trees are ultrametric
    #tree<-sim.bd.taxa.age(ntax,nsim,lambda,mu,frac,age=age)
    tree<-sim.bd.taxa(n=ntax,numbsim=nsim,lambda=lambda,mu=mu,complete = FALSE, stochsampling = TRUE)
    
    #tree<-sim.bd.taxa(n=ntax,nsim,lambda,mu,frac)
    #plot(tree[[1]])
    sim.taxa.kappa.array[sizeIndex,ageIndex,1]<-ntax
    vcvtree<-vcv(tree[[1]])
    vcvtree<-vcvtree
    sim.taxa.kappa.array[sizeIndex,ageIndex,2]<-kappa(vcvtree/max(vcvtree))
  }
}

dim(sim.taxa.kappa.array)
fullsim.taxa.kappa.array<-sim.taxa.kappa.array
sim.taxa.kappa.array<-apply(sim.taxa.kappa.array,c(1,3),mean)
sim.taxa.kappa.array<- data.frame(taxa=sim.taxa.kappa.array[,1],kappa=sim.taxa.kappa.array[,2])

# 4. Take into account tree shape.
# 4.1 Are more balanced trees more stable than unbalanced tips?

# 4.2 What about _tipness_? Are more tipy trees less stable ?  
#     trees that show a more distributed brnaching pattern from root to the tips? 

lm_rawtree<- lm(log10(kappa) ~ taxa, data=taxa.kappa.array)
lm_simtree<- lm(log10(kappa) ~ taxa, data=sim.taxa.kappa.array)

print(lm_rawtree)
print(lm_simtree)


library(gridExtra)
#p1<-ggplot(data=taxa.kappa.array,mapping = aes(x=taxa,y=log10(kappa))) + geom_point( color = "blue") + 
#  geom_smooth(method=lm, color="black") +
#  labs(title="          kappa vs. taxa (Real Data)", y = expression(paste(log[10], " ", kappa)), x="number of taxa")  + theme_classic()
#p1<-p1 + scale_y_log10(limits=c(2, 10))   #+scale_y_continuous(trans = 'log10',limits=c(3, 17))


p1<- ggplot(data=taxa.kappa.array,mapping = aes(x=taxa,y=log10(kappa))) + geom_point( color = "blue") + 
  geom_smooth(method=lm, se=TRUE, color="black") +
  labs(title="          kappa vs. taxa (Real Data)", y = expression(paste(log[10], " ", kappa)), x="number of taxa")  + theme_classic()

p1<-p1 + scale_y_log10(limits=c(1, 7))

p2<-ggplot(data=sim.taxa.kappa.array,mapping = aes(x=taxa,y=log10(kappa))) + geom_point( color = "blue") + 
  geom_smooth(method=lm, se=TRUE, color="black") +
  labs(title="          kappa vs. taxa (Real Data)", y = expression(paste(log[10], " ", kappa)), x="number of taxa")  + theme_classic()

p2<-p2 + scale_y_log10(limits=c(1, 7))


grid.arrange(p1, p2, nrow = 1)


#save.image("fig2.rda")
