rm(list=ls())
setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaShrinkSummary/")

library(ggplot2)
library(devtools)
library(qpcR)
#install_github("kassambara/easyGgplot2")
library(easyGgplot2)
require(gridExtra)

foldername<-paste("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaShrinkRun",c(1,2,3),sep="")
mu.array=c(0,10)
sigmasq.array=c(1,5)
treesize.array<-c(100,500)
delta.array<-c(0.001,0.005,0.01,1)

for(paraIndex in 1:2){   
#  paraIndex<-1
  mu<-mu.array[paraIndex]
  sigmasq<-sigmasq.array[paraIndex]
  for(sizeIndex in seq_along(treesize.array)){
#    sizeIndex<-1
    treesize<-treesize.array[sizeIndex]
    solfiletoload<-paste("solshrinksize",treesize,"mu",mu,"sigmasq",sigmasq,".rda",sep="")
    solve_shrink_raw_delta_combined<-NULL
    solve_shrink_RMSD_delta_0.001_combined<-NULL
    solve_shrink_RMSD_delta_0.005_combined<-NULL
    solve_shrink_RMSD_delta_0.01_combined<-NULL
    solve_shrink_RMSD_delta_1_combined<-NULL
    
    unsolfiletoload<-paste("unsolshrinksize",treesize,"mu",mu,"sigmasq",sigmasq,".rda",sep="")
    unsolve_shrink_RMSD_delta_0.001_combined<-NULL
    unsolve_shrink_RMSD_delta_0.005_combined<-NULL
    unsolve_shrink_RMSD_delta_0.01_combined<-NULL
    unsolve_shrink_RMSD_delta_1_combined<-NULL
    
    for(folderIndex in seq_along(foldername)){
#      folderIndex<-1
      filepath<-foldername[folderIndex]
      setwd(filepath)
      load(solfiletoload)
      solve_shrink_raw_delta_combined<-rbind(solve_shrink_raw_delta_combined, solve_shrink_raw_delta0.001,solve_shrink_raw_delta0.005,solve_shrink_raw_delta0.01,solve_shrink_raw_delta1)
      solve_shrink_RMSD_delta_0.001_combined<-rbind(solve_shrink_RMSD_delta_0.001_combined, solve_shrink_RMSD_delta0.001)
      solve_shrink_RMSD_delta_0.005_combined<-rbind(solve_shrink_RMSD_delta_0.005_combined, solve_shrink_RMSD_delta0.005)
      solve_shrink_RMSD_delta_0.01_combined<-rbind(solve_shrink_RMSD_delta_0.01_combined, solve_shrink_RMSD_delta0.01)
      solve_shrink_RMSD_delta_1_combined<-rbind(solve_shrink_RMSD_delta_1_combined, solve_shrink_RMSD_delta1)
    
      load(unsolfiletoload)
      unsolve_shrink_RMSD_delta_0.001_combined<-rbind(unsolve_shrink_RMSD_delta_0.001_combined, unsolve_shrink_RMSD_delta0.001)
      unsolve_shrink_RMSD_delta_0.005_combined<-rbind(unsolve_shrink_RMSD_delta_0.005_combined, unsolve_shrink_RMSD_delta0.005)
      unsolve_shrink_RMSD_delta_0.01_combined<-rbind(unsolve_shrink_RMSD_delta_0.01_combined, unsolve_shrink_RMSD_delta0.01)
      unsolve_shrink_RMSD_delta_1_combined<-rbind(unsolve_shrink_RMSD_delta_1_combined, unsolve_shrink_RMSD_delta1)
      }

    
    dimraw<-dim(solve_shrink_raw_delta_combined)[1]
    dim001<-dim(solve_shrink_RMSD_delta_0.001_combined)[1]
    dim005<-dim(solve_shrink_RMSD_delta_0.005_combined)[1]
    dim01<-dim(solve_shrink_RMSD_delta_0.01_combined)[1]
    dim1<-dim(solve_shrink_RMSD_delta_1_combined)[1]
    df.plot.sol<-data.frame(RMSD=c(solve_shrink_raw_delta_combined[,2],
                                 solve_shrink_RMSD_delta_0.001_combined[,2],
                                 solve_shrink_RMSD_delta_0.005_combined[,2],
                                 solve_shrink_RMSD_delta_0.01_combined[,2],
                                 solve_shrink_RMSD_delta_1_combined[,2]),
                            Shrinkage = c(
                            rep("raw",dimraw),
                            rep("δ=0.001",dim001),
                            rep("δ=0.005",dim005),
                            rep("δ=0.01",dim01),
                            rep("δ=1",dim1)
                          ),
                          Taxa = c(
                            rep(treesize, dimraw+dim001+dim005+dim01+dim1)
                          ))
    setwd("~/Dropbox/FCU/Teaching/Mentoring/2019Spring/ChienSanPeng/KappaShrinkSummary/")
    
    assign(paste( "dfplot.sol",treesize ,sep=""), df.plot.sol)
    
    dim001<-dim(unsolve_shrink_RMSD_delta_0.001_combined)[1]
    dim005<-dim(unsolve_shrink_RMSD_delta_0.005_combined)[1]
    dim01<-dim(unsolve_shrink_RMSD_delta_0.01_combined)[1]
    dim1<-dim(unsolve_shrink_RMSD_delta_1_combined)[1]
    df.plot.unsol<-data.frame(RMSD=c(
                                   unsolve_shrink_RMSD_delta_0.001_combined[,2],
                                   unsolve_shrink_RMSD_delta_0.005_combined[,2],
                                   unsolve_shrink_RMSD_delta_0.01_combined[,2],
                                   unsolve_shrink_RMSD_delta_1_combined[,2]),
                            Shrinkage = c(
                              rep("δ=0.001",dim001),
                              rep("δ=0.005",dim005),
                              rep("δ=0.01",dim01),
                              rep("δ=1",dim1)
                            ),
                            Taxa = c(
                              rep(treesize, dimraw+dim001+dim005+dim01+dim1)
                            ))
    
    assign(paste( "dfplot.unsol",treesize ,sep=""), df.plot.unsol)
    
    }#sizeIndex
   
  df.solsigmasq<-rbind(dfplot.sol100,dfplot.sol500)
  assign(paste("df.solsigmasq",sigmasq,sep=""),df.solsigmasq)
  df.solsigmasq$Taxa <- factor(df.solsigmasq$Taxa, levels = c("100","500"))
  df.solsigmasq$Shrinkage<-factor(df.solsigmasq$Shrinkage, levels = c("raw", "δ=0.001","δ=0.005","δ=0.01","δ=1")) 
  p.sol<-ggplot(df.solsigmasq,aes(x=Taxa,y=RMSD)) + geom_boxplot(aes(fill=Shrinkage)) + scale_y_continuous(trans = "log10")
  p.sol<- p.sol + ggtitle(label = paste("σ = ", sigmasq, ", Sol  C",sep=""))#, subtitle = paste("σ2 = ", sigmasq , ", n = " , size,sep=""))
  p.sol<- p.sol + theme(plot.title = element_text(color = "black", size = 20, face = "bold")) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste("p.sol.sigmasq",sigmasq,sep=""),p.sol)
  
  df.unsolsigmasq<-rbind(dfplot.unsol100,dfplot.unsol500)
  assign(paste("df.unsolsigmasq",sigmasq,sep=""),df.unsolsigmasq)
  df.unsolsigmasq$Taxa <- factor(df.unsolsigmasq$Taxa, levels = c("100","500"))
  df.unsolsigmasq$Shrinkage<-factor(df.unsolsigmasq$Shrinkage, levels = c("raw", "δ=0.001","δ=0.005","δ=0.01","δ=1")) 
  p.unsol<-ggplot(df.unsolsigmasq,aes(x=Taxa,y=RMSD)) + geom_boxplot(aes(fill=Shrinkage)) + scale_y_continuous(trans = "log10")
  p.unsol<- p.unsol + ggtitle(label = paste("σ = ", sigmasq, ", Unsol C", sep=""))#, subtitle = paste("σ2 = ", sigmasq , ", n = " , size,sep=""))
  p.unsol<- p.unsol + theme(plot.title = element_text(color = "black", size = 20, face = "bold")) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste("p.unsol.sigmasq",sigmasq,sep=""),p.unsol)
  }#paraIndex

#Sys.sleep(time=8)
grid.arrange(p.sol.sigmasq1,p.sol.sigmasq5,p.unsol.sigmasq1,p.unsol.sigmasq5,ncol=2)


    


