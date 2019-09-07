setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaShrinkRun3/")
rm(list=ls())
treesize.array<-c(100, 500)#, 800)
mu.array<-c(0,10) # or mu <- 10
sigmasq.array<-c(1,5)# sigmasq<- 5
cases<-c("sol","unsol")

for(treesizeIndex in 1:length(treesize.array)){
  treesize<-treesize.array[treesizeIndex]
  for(caseIndex in 1:2){
    case<-cases[caseIndex]
  for(paraIndex in 1:2){
    mu<-mu.array[paraIndex]
    sigmasq<-sigmasq.array[paraIndex]
    jobname<-paste(case,"shrinksize",treesize,"mu",mu,"sigmasq",sigmasq,".r",sep="")
    print(jobname)
    jobgo<-paste("nohup R CMD BATCH ", jobname , " >/dev/null &", sep="") 
    print(jobgo)
    system(jobgo)
        }
  }
}

