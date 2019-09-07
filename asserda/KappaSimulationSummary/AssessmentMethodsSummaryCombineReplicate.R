rm(list=ls())
setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationSummary")
library(ggplot2)
library(devtools)
library(qpcR)
#install_github("kassambara/easyGgplot2")
library(easyGgplot2)
require(gridExtra)

foldername<-paste("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationRun",c(1,2,3,4,5,6),sep="")

mu.array=c(0,10)
sigmasq.array=c(1,5)
treesize.array<-c(100,500)

for(sizeIndex in seq_along(treesize.array)){
   size<-treesize.array[sizeIndex]
   for(paraIndex in 1:2){   
       print(paste("size", size, "paraset",paraIndex,sep=""))
       mu<-mu.array[paraIndex]
       sigmasq<-sigmasq.array[paraIndex]
       
   filetoload<-paste("size",size,"mu",mu, "sigmasq",sigmasq,".rda",sep="")

solve_raw_RMSD_combined<-NULL
solve_shrink_RMSD_combined<-NULL
solve_stretch_RMSD_combined<-NULL
solve_droptip_RMSD_combined<-NULL
solve_prun_RMSD_combined<-NULL

unsolve_shrink_RMSD_combined<-NULL
unsolve_stretch_RMSD_combined<-NULL
unsolve_droptip_RMSD_combined<-NULL
unsolve_prun_RMSD_combined<-NULL

for(folderIndex in seq_along(foldername)){
    solve_raw_RMSD<-NULL
    solve_shrink_RMSD<-NULL
    solve_stretch_RMSD<-NULL
    solve_droptip_RMSD<-NULL
    solve_prun_RMSD<-NULL
    
    unsolve_shrink_RMSD<-NULL
    unsolve_stretch_RMSD<-NULL
    unsolve_droptip_RMSD<-NULL
    unsolve_prun_RMSD<-NULL
#    folderIndex<-1
    folderpath<-foldername[folderIndex]
    setwd(folderpath)  
    try(load(filetoload))
    try(solve_raw_RMSD_combined<-rbind(solve_raw_RMSD_combined,solve_raw_RMSD))
    try(solve_shrink_RMSD_combined<-rbind(solve_shrink_RMSD_combined,solve_shrink_RMSD))
    try(solve_stretch_RMSD_combined<-rbind(solve_stretch_RMSD_combined,solve_stretch_RMSD))
    try(solve_droptip_RMSD_combined<-rbind(solve_droptip_RMSD_combined,solve_droptip_RMSD))
    try(solve_prun_RMSD_combined<-rbind(solve_prun_RMSD_combined,solve_prun_RMSD))
    
    try(unsolve_shrink_RMSD_combined<-rbind(unsolve_shrink_RMSD_combined,unsolve_shrink_RMSD))
    try(unsolve_stretch_RMSD_combined<-rbind(unsolve_stretch_RMSD_combined,unsolve_stretch_RMSD))
    try(unsolve_droptip_RMSD_combined<-rbind(unsolve_droptip_RMSD_combined,unsolve_droptip_RMSD))
    try(unsolve_prun_RMSD_combined<-rbind(unsolve_prun_RMSD_combined,unsolve_prun_RMSD))
    }

savefilename<-paste("size",size,"mu",mu, "sigmasq",sigmasq,"Combined.rda",sep="")
setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaSimulationSummary")
save.image(savefilename)


# To remove outlier for better boxplot
# solve_raw_RMSD_combined<-TheCutoffSet(solve_raw_RMSD_combined,L=L,U=U)
# solve_shrink_RMSD_combined<-TheCutoffSet(solve_shrink_RMSD_combined,L=L,U=U)
# solve_stretch_RMSD_combined<-TheCutoffSet(solve_stretch_RMSD_combined,L=L,U=U)
# solve_droptip_RMSD_combined<-TheCutoffSet(solve_droptip_RMSD_combined,L=L,U=U)    
# solve_prun_RMSD_combined<-TheCutoffSet(solve_prun_RMSD_combined,L=L,U=U)
# 
# unsolve_shrink_RMSD_combined<-TheCutoffSet(unsolve_shrink_RMSD_combined,L=L,U=U)
# unsolve_stretch_RMSD_combined<-TheCutoffSet(unsolve_stretch_RMSD_combined,L=L,U=U)
# unsolve_droptip_RMSD_combined<-TheCutoffSet(unsolve_droptip_RMSD_combined,L=L,U=U)
# unsolve_prun_RMSD_combined<-TheCutoffSet(unsolve_prun_RMSD_combined,L=L,U=U)

try(solve_raw_RMSD_combined<-solve_raw_RMSD_combined[!is.nan(solve_raw_RMSD_combined[,1]),])
try(solve_shrink_RMSD_combined<-solve_shrink_RMSD_combined[!is.nan(solve_shrink_RMSD_combined[,1]),])
try(solve_stretch_RMSD_combined<-solve_stretch_RMSD_combined[!is.nan(solve_stretch_RMSD_combined[,1]),])
try(solve_droptip_RMSD_combined<-solve_droptip_RMSD_combined[!is.nan(solve_droptip_RMSD_combined[,1]),])
try(solve_prun_RMSD_combined<-solve_prun_RMSD_combined[!is.nan(solve_prun_RMSD_combined[,1]),])

try(unsolve_shrink_RMSD_combined<-unsolve_shrink_RMSD_combined[!is.nan(unsolve_shrink_RMSD_combined[,1]),])
try(unsolve_stretch_RMSD_combined<-unsolve_stretch_RMSD_combined[!is.nan(unsolve_stretch_RMSD_combined[,1]),])
try(unsolve_droptip_RMSD_combined<-unsolve_droptip_RMSD_combined[!is.nan(unsolve_droptip_RMSD_combined[,1]),])
try(unsolve_prun_RMSD_combined<-unsolve_prun_RMSD_combined[!is.nan(unsolve_prun_RMSD_combined[,1]),])


###MAKE DATAFRAME

solraw<-data.frame(solraw=solve_raw_RMSD_combined[,1])
solshrink<-data.frame(solshrink=solve_shrink_RMSD_combined[,1])
solstretch<-data.frame(solstretch=solve_stretch_RMSD_combined[,1])
soldroptip<-data.frame(soldroptip=solve_droptip_RMSD_combined[,1])
solprun<-data.frame(solprun=solve_prun_RMSD_combined[,1])
unsolshrink<-data.frame(unsolshrink=unsolve_shrink_RMSD_combined[,1])
unsolstretch<-data.frame(unsolstretch=unsolve_stretch_RMSD_combined[,1])
unsoldroptip<-data.frame(unsoldroptip=unsolve_droptip_RMSD_combined[,1])
unsolprun<-data.frame(unsolprun=unsolve_prun_RMSD_combined[,1])

df.theta<-qpcR:::cbind.na(solraw,solshrink,solstretch,soldroptip,solprun,unsolshrink,unsolstretch,unsoldroptip,unsolprun)

head(df.theta)
apply(df.theta,2,summary)

df.theta.plot <-  data.frame(RMSD=c(df.theta$solraw,df.theta$solshrink,df.theta$solstretch,df.theta$soldroptip,df.theta$solprun,  rep(0,length(df.theta$solraw)), df.theta$unsolshrink,df.theta$unsolstretch,df.theta$unsoldroptip,df.theta$unsolprun),
                             Method= c(  rep("sol", length(df.theta$solraw)),  
                                         rep("sol", length(df.theta$solshrink)),  
                                         rep("sol", length(df.theta$solstretch)),  
                                         rep("sol", length(df.theta$soldroptip)),  
                                         rep("sol", length(df.theta$solprun)),  

                                         rep("unsol", length(df.theta$solraw)),
                                         rep("unsol", length(df.theta$unsolshrink)),  
                                         rep("unsol", length(df.theta$unsolstretch)),  
                                         rep("unsol", length(df.theta$unsoldroptip)),
                                         rep("unsol", length(df.theta$unsolprun))),
                             Type=c(     
                                         rep("raw", length(df.theta$solraw)),
                                         rep("shrink", length(df.theta$solshrink)),
                                         rep("stretch", length(df.theta$solstretch)), 
                                         rep("droptip", length(df.theta$soldroptip)),
                                         rep("prun", length(df.theta$unsolprun)),
                                         rep("raw", length(df.theta$solraw)),
                                         rep("shrink", length(df.theta$unsolshrink)),
                                         rep("stretch", length(df.theta$unsolstretch)),
                                         rep("droptip", length(df.theta$unsoldroptip)),
                                         rep("prun", length(df.theta$solprun))
                             ))     

df.theta.plot$Type <- factor(df.theta.plot$Type, levels = c("droptip", "prun", "raw","stretch","shrink"))


### MAKE PLOT
#df.theta.plot<-df.theta.plot[order(df.theta.plot$RMSD),]
#df.theta.plot$Method <- factor(df.theta.plot$Method,levels=df.theta.plot$Method[order(df.theta.plot$RMSD)])

head(df.theta.plot)
tail(df.theta.plot)

#p.theta<-ggplot(df.theta.plot,aes(x=Type,y=RMSD)) + geom_boxplot(aes(fill=Method)) +scale_y_continuous(trans = "log10")

#p.theta

#p.violin.theta<-ggplot2.violinplot(data=df.theta.plot,xName="method",yName="rmsd")
#p.violin.theta

#### sigma sq

solraw<-data.frame(solraw=solve_raw_RMSD_combined[,2])
solshrink<-data.frame(solshrink=solve_shrink_RMSD_combined[,2])
solstretch<-data.frame(solstretch=solve_stretch_RMSD_combined[,2])
soldroptip<-data.frame(soldroptip=solve_droptip_RMSD_combined[,2])
solprun<-data.frame(solprun=solve_prun_RMSD_combined[,2])
unsolshrink<-data.frame(unsolshrink=unsolve_shrink_RMSD_combined[,2])
unsolstretch<-data.frame(unsolstretch=unsolve_stretch_RMSD_combined[,2])
unsoldroptip<-data.frame(unsoldroptip=unsolve_droptip_RMSD_combined[,2])
unsolprun<-data.frame(unsolprun=unsolve_prun_RMSD_combined[,2])


df.sigmasq<-qpcR:::cbind.na(solraw,solshrink,solstretch,soldroptip,solprun,unsolshrink,unsolstretch,unsoldroptip,unsolprun)

head(df.sigmasq)
apply(df.sigmasq,2,summary)

df.sigmasq.plot <-  data.frame(RMSD=c(df.sigmasq$solraw,df.sigmasq$solshrink,df.sigmasq$solstretch,df.sigmasq$soldroptip,df.sigmasq$solprun,  rep(0,length(df.sigmasq$solraw)), df.sigmasq$unsolshrink,df.sigmasq$unsolstretch,df.sigmasq$unsoldroptip,df.sigmasq$unsolprun),
                               Method= c(  rep("sol", length(df.sigmasq$solraw)),  
                                           rep("sol", length(df.sigmasq$solshrink)),  
                                           rep("sol", length(df.sigmasq$solstretch)),  
                                           rep("sol", length(df.sigmasq$soldroptip)),  
                                           rep("sol", length(df.sigmasq$solprun)),  
                                           
                                           rep("unsol", length(df.sigmasq$solraw)),
                                           rep("unsol", length(df.sigmasq$unsolshrink)),  
                                           rep("unsol", length(df.sigmasq$unsolstretch)),  
                                           rep("unsol", length(df.sigmasq$unsoldroptip)),
                                           rep("unsol", length(df.sigmasq$unsolprun))),
                               Type=c(     
                                   rep("raw", length(df.sigmasq$solraw)),
                                   rep("shrink", length(df.sigmasq$solshrink)),
                                   rep("stretch", length(df.sigmasq$solstretch)), 
                                   rep("droptip", length(df.sigmasq$soldroptip)),
                                   rep("prun", length(df.sigmasq$unsolprun)),
                                   rep("raw", length(df.sigmasq$solraw)),
                                   rep("shrink", length(df.sigmasq$unsolshrink)),
                                   rep("stretch", length(df.sigmasq$unsolstretch)),
                                   rep("droptip", length(df.sigmasq$unsoldroptip)),
                                   rep("prun", length(df.sigmasq$solprun))
                               ))     
df.sigmasq.plot$Type <- factor(df.sigmasq.plot$Type, levels = c("droptip", "prun", "raw","stretch","shrink"))

### MAKE PLOT
#df.sigmasq.plot<-df.sigmasq.plot[order(df.sigmasq.plot$RMSD),]
#df.sigmasq.plot$Method <- factor(df.sigmasq.plot$Method,levels=df.sigmasq.plot$Method[order(df.sigmasq.plot$RMSD)])

head(df.sigmasq.plot)
tail(df.sigmasq.plot)

p.sigmasq<-ggplot(df.sigmasq.plot,aes(x=Type,y=RMSD)) + geom_boxplot(aes(fill=Method))  +scale_y_continuous(trans = "log10") 
p.sigmasq<- p.sigmasq + ggtitle(label = paste("σ2 = ", sigmasq , ", n = " , size,sep=""))#, subtitle = paste("σ2 = ", sigmasq , ", n = " , size,sep=""))
p.sigmasq<- p.sigmasq + theme(plot.title = element_text(color = "black", size = 20, face = "bold")) + theme(plot.title = element_text(hjust = 0.5))
#    plot.subtitle = element_text(color = "blue")


p.theta<-ggplot(df.theta.plot,aes(x=Type,y=RMSD)) + geom_boxplot(aes(fill=Method)) +scale_y_continuous(trans = "log10")
p.theta<- p.theta + ggtitle(label = paste("θ = ", mu , ", n = ", size ,sep=""))#, subtitle = paste("θ = ", mu , ", n = ", size ,sep=""))
p.theta<- p.theta + theme(
    plot.title = element_text(color = "black", size = 20, face = "bold"))+ theme(plot.title = element_text(hjust = 0.5))
#    plot.subtitle = element_text(color = "blue")

assign(paste("pt.theta",mu,"sigmasq",sigmasq,"size",size,sep=""), p.theta)
assign(paste("ps.theta",mu,"sigmasq",sigmasq,"size",size,sep=""), p.sigmasq)

 }
}#p.sigmasq


grid.arrange(pt.theta0sigmasq1size100,pt.theta10sigmasq5size100,pt.theta0sigmasq1size500,pt.theta10sigmasq5size500,ncol=2)
grid.arrange(ps.theta0sigmasq1size100,ps.theta10sigmasq5size100,ps.theta0sigmasq1size500,ps.theta10sigmasq5size500,ncol=2)

#grid.arrange(pt.theta0sigmasq1size100,pt.theta10sigmasq5size100,pt.theta0sigmasq1size500,pt.theta10sigmasq5size500,pt.theta10sigmasq5size800,ncol=3)
#grid.arrange(ps.theta0sigmasq1size100,ps.theta10sigmasq5size100,ps.theta0sigmasq1size500,ps.theta10sigmasq5size500,ps.theta10sigmasq5size800,ncol=3)

#grid.arrange(pt.theta0sigmasq1size500,pt.theta10sigmasq5size500,ncol=2)
#grid.arrange(ps.theta0sigmasq1size100,ps.theta10sigmasq5size100,ncol=2)


