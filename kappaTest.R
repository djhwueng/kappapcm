
setwd("~/Dropbox/CollabAdamsCollyerJhwuengOMeara/")
library(TreeSim)
library(ape)
library(MuMIn)
install.packages("MuMIn")
nsim = 100
nreps = 200


GetSummary <- function(phy, lambda=NA, mu=NA, frac=NA, age=NA) {
	log.kappa.result <- log(kappa(vcv(phy, exact=TRUE)))
	min.brlen <- min(phy$edge.length)
	max.brlen <- max(phy$edge.length)
	brlen.ratio <- max.brlen / min.brlen
	brlen.var <- var(phy$edge.length)
	brlen.median <- median(phy$edge.length)
	tip.lengths <- phy$edge.length[which(phy$edge[,2]<=Ntip(phy))]
	internal.lengths <- phy$edge.length[which(phy$edge[,2]>Ntip(phy))]
	min.tip <- min(tip.lengths)
	max.tip <- max(tip.lengths)
	min.internal <- min(internal.lengths)
	max.internal <- max(internal.lengths)
	return(data.frame(Ntip=Ntip(phy), log.kappa=log.kappa.result, min.brlen=min.brlen, max.brlen=max.brlen, brlen.ratio=brlen.ratio, brlen.var=brlen.var, brlen.median=brlen.median, min.tip=min.tip, max.tip=max.tip, min.internal=min.internal, max.internal=max.internal, lambda=lambda, mu=mu, frac=frac, age=age))
}


all.results <- data.frame()
for (i in sequence(nreps)) {
	ntax = round(runif(1, 10, 100))
	lambda = runif(1, 0.01, 0.1)
	mu= runif(1, 0, lambda)
	frac= runif(1, 0.1, 1)
	age=log(ntax) / (lambda-mu)
	trees <- sim.bd.taxa.age(ntax, nsim, lambda, mu, frac, age)
	tree.results <- t(sapply(trees, GetSummary, lambda=lambda, mu=mu, frac=frac, age=age))
	if (i==1) {
		all.results <- tree.results
	} else {
		all.results <- rbind(all.results, tree.results)
	}
	print(paste("Done",i,"of",nreps))
}
all.results <- data.frame(all.results)
for (i in sequence(dim(all.results)[2])) {
	all.results[,i] <- unlist(all.results[,i])	
}

#save(list=ls(), file="~/Desktop/kappaTest.rda")


load("kappaTest.rda")

options(na.action = "na.fail")
global.model <- lm(log.kappa ~ ., data=all.results)
dredge.results <- dredge(global.model, m.lim=c(NA,4), trace=1)
print(dredge.results, abbrev.names=FALSE)
head(dredge.results)
write.csv(dredge.results, file="~/Desktop/kappaTest.csv")
save(list=ls(), file="~/Desktop/kappaTest.rda")