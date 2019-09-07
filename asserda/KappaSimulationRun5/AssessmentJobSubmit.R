setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationRun5/")
rm(list=ls())
treesize.array<-c(100, 500, 800)
mu.array<-c(0,10) # or mu <- 10
sigmasq.array<-c(1,5)# sigmasq<- 5

for(treesizeIndex in 1:length(treesize.array)){
  treesize<-treesize.array[treesizeIndex]
  for(paraIndex in 1:2){
    mu<-mu.array[paraIndex]
    sigmasq<-sigmasq.array[paraIndex]
    jobname<-paste("size",treesize,"mu",mu,"sigmasq",sigmasq,".r",sep="")
    print(jobname)
    jobgo<-paste("nohup R CMD BATCH ", jobname , " >/dev/null &", sep="") 
    print(jobgo)
    system(jobgo)
        }
  }



