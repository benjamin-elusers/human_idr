source("/Users/atarg/Documents/Levy_lab/scripts/RToolBox2/MicroscopeToolBox.R")
options(scipen=999)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(car)
library(effsize)
library(pwr)
library(dplyr)
library(pheatmap)
library(Rfast)
library(forcats)
library(gridExtra)
library(grid)
library(reshape2)
library(igraph)
library(RColorBrewer)
library(viridis)
library(ggraph)
library(umap)
library(ggfortify)
library(cluster)

#Screen directory:
screen.dir=
  "/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/"


#Read .RDS raw data results files:
output.dir=
  paste0(screen.dir, "results/")

data.KPlib=
  readRDS(paste0(output.dir,"raw_list.RDS"))[[1]]
design.KPlib=
  readRDS(paste0(output.dir,"design.RDS"))

#function to check which pictures (s_###) belong to a well in screen:
# get.sNUM(design.KPlib,plate=1,well="M14")



# Call on functions to correct raw data:

data.KPlib=filter.saturated(data.KPlib,"GFP_int_b0")
data.KPlib=filter.saturated(data.KPlib,"RFP_int_b0")

data.KPlib=correct.function(data.KPlib,"GFP_int_b5","RFP_int_b5")
data.KPlib=correct.function(data.KPlib,"GFP_int_b0","RFP_int_b0")
data.KPlib=correct.function(data.KPlib,"GFP_int_mean","RFP_int_mean")
data.KPlib=correct.function(data.KPlib,"f1_inGFP_toGFPmean","f1_inGFP_toRFPmean")

platedef=
  read.csv(paste0(screen.dir,"platedef_idr.csv"))
colnames(platedef)=c("WELL","BF","GFP","RFP")
GFP.list.full=c("CycT1","Ddx3x","FUS","GATA3","G3BP1","hnRNPA1",
                #"HP1a",
                "HspB8","Nephrin","NPM1","RBM14","Syn1",
                "TAF15","TDP43","Valency","Syn2")
RFP.list.full=c("CycT1","Ddx3x","FUS","GATA3","G3BP1","hnRNPA1",
                "HP1a","HspB8","Nephrin","NPM1","RBM14","Syn1",
                "TAF15","TDP43","Valency","Ddx4","DYRK3","Era","NUP98",
                "NUP98_trunc","UBQ2","Syn2")




############################################################################

GFP.valency.plot.list <- list()
GFP.csat.plot.list <- list()
GFP.frac.plot.list <- list()
sat.val.GFP=c()
sat.qual.GFP=c()

for (I in c(1:length(GFP.list.full))){
  G = GFP.list.full[I]
  R = "Valency"
  W = platedef$WELL[which(platedef$GFP==G & platedef$RFP == R)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("YFP-",G)
  dtmp$GFP_RFP_b5=as.numeric(dtmp$RFP_int_mean)
  if (ncell>0){
    
    q1=quantile(dtmp$GFP_RFP_b5,prob=0.20)
    q2=quantile(dtmp$GFP_RFP_b5,prob=0.40)
    q3=quantile(dtmp$GFP_RFP_b5,prob=0.60)
    q4=quantile(dtmp$GFP_RFP_b5,prob=0.80)
    punc.GFP= (dtmp$inGFPnfoci>0 & dtmp$f1_inGFP_toGFPmean>150)
    
    data.q1= dtmp %>% filter(GFP_RFP_b5<=q1 & punc.GFP)
    data.q2= dtmp %>% filter(GFP_RFP_b5<=q2 & GFP_RFP_b5>q1 & punc.GFP)
    data.q3= dtmp %>% filter(GFP_RFP_b5<=q3 & GFP_RFP_b5>q2 & punc.GFP)
    data.q4= dtmp %>% filter(GFP_RFP_b5<=q4 & GFP_RFP_b5>q3 & punc.GFP)
    data.q5= dtmp %>% filter(GFP_RFP_b5>q4 & punc.GFP)
    
    data=list(data.q1,data.q2,data.q3,data.q4,data.q5)
    
    frac.punc.q1 = nrow(data.q1)/sum(dtmp$GFP_RFP_b5<q1)
    frac.punc.q2 = nrow(data.q2)/sum(dtmp$GFP_RFP_b5<q2 & dtmp$GFP_RFP_b5>q1)
    frac.punc.q3 = nrow(data.q3)/sum(dtmp$GFP_RFP_b5<q3 & dtmp$GFP_RFP_b5>q2)
    frac.punc.q4 = nrow(data.q4)/sum(dtmp$GFP_RFP_b5<q4 & dtmp$GFP_RFP_b5>q3)
    frac.punc.q5 = nrow(data.q5)/sum(dtmp$GFP_RFP_b5>q4)
    #round the frac.punc list!!!##### round(frac.punc[I],3)
    frac.punc.list=list(round(frac.punc.q1,3),round(frac.punc.q2,3),round(frac.punc.q3,3),round(frac.punc.q4,3),round(frac.punc.q5,3))
    
    if (nrow(data.q1)>20 & nrow(data.q2)>20 & nrow(data.q3)>20 & nrow(data.q4)>20 & nrow(data.q5)>20){
      SUM=fig.data(n=5,FP="GFP")
      val.plot=val.plot.full(n=5,FP="GFP")
      GFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=5)
      csat=csat.plot(FP="GFP")
      GFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=5)
      frac.plot=frac.plot.full(FP="GFP")
      GFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20 & nrow(data.q2)>20 & nrow(data.q3)>20 & nrow(data.q4)>20){
      SUM=fig.data(n=4,FP="GFP")
      val.plot=val.plot.full(n=4,FP="GFP")
      GFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=4)
      csat=csat.plot(FP="GFP")
      GFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=4)
      frac.plot=frac.plot.full(FP="GFP")
      GFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20 & nrow(data.q2)>20 & nrow(data.q3)>20){
      SUM=fig.data(n=3,FP="GFP")
      val.plot=val.plot.full(n=3,FP="GFP")
      GFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=3)
      csat=csat.plot(FP="GFP")
      GFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=3)
      frac.plot=frac.plot.full(FP="GFP")
      GFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20 & nrow(data.q2)>20) {
      SUM=fig.data(n=2,FP="GFP")
      val.plot=val.plot.full(n=2,FP="GFP")
      GFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=2)
      csat=csat.plot(FP="GFP")
      GFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=2)
      frac.plot=frac.plot.full(FP="GFP")
      GFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20) {
      SUM=fig.data(n=1,FP="GFP")
      val.plot=val.plot.full(n=1,FP="GFP")
      GFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=1)
      csat=csat.plot(FP="GFP")
      GFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=1)
      frac.plot=frac.plot.full(FP="GFP")
      GFP.frac.plot.list[[I]]=frac.plot
      
    } else {
      val.plot=val.plot.full(n=0,FP="GFP")
      GFP.valency.plot.list[[I]]=val.plot

      csat.data.df=csat.data(n=0)
      csat=csat.plot(FP="GFP")
      GFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=0)
      frac.plot=frac.plot.full(FP="GFP")
      GFP.frac.plot.list[[I]]=frac.plot
    }
  }
  else{
    val.plot=val.plot.full(n=0,FP="GFP")
    GFP.valency.plot.list[[I]]=val.plot
    
    csat.data.df=csat.data(n=0)
    csat=csat.plot(FP="GFP")
    GFP.csat.plot.list[[I]]=csat
    
    frac.data.df=frac.data(n=0)
    frac.plot=frac.plot.full(FP="GFP")
    GFP.frac.plot.list[[I]]=frac.plot
    
  }
  #print(GFP.valency.plot.list[[I]])
  #print(GFP.csat.plot.list[[I]])
  #print(GFP.frac.plot.list[[I]])
}


panel.GFP=ggarrange(plotlist=GFP.valency.plot.list, ncol=1,nrow=16)
panel.sat=ggarrange(plotlist=GFP.csat.plot.list, ncol=1,nrow=16)
panel.frac=ggarrange(plotlist=GFP.frac.plot.list, ncol=1,nrow=16)

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/val_2_new.pdf",width=8,height=128)
print(panel.GFP)
dev.off()

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/csat_2_new.pdf",width=8,height=128)
print(panel.sat)
dev.off()

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/frac_2_new.pdf",width=8,height=128)
print(panel.frac)
dev.off()
######################################################


RFP.valency.plot.list <- list()
RFP.csat.plot.list <- list()
RFP.frac.plot.list <- list()
sat.val.RFP=c()
sat.qual.RFP=c()

for (I in c(1:length(RFP.list.full))){
  G = "Valency"
  R = RFP.list.full[I]
  W = platedef$WELL[which(platedef$GFP==G & platedef$RFP == R)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("RFP-",R)
  dtmp$GFP_RFP_b5=as.numeric(dtmp$GFP_int_mean)
  if (ncell>0){
    
    q1=quantile(dtmp$GFP_RFP_b5,prob=0.20)
    q2=quantile(dtmp$GFP_RFP_b5,prob=0.40)
    q3=quantile(dtmp$GFP_RFP_b5,prob=0.60)
    q4=quantile(dtmp$GFP_RFP_b5,prob=0.80)
    punc.RFP= (dtmp$inRFPnfoci>0 & dtmp$f1_inRFP_toRFPmean>150)
    
    data.q1= dtmp %>% filter(GFP_RFP_b5<=q1 & punc.RFP)
    data.q2= dtmp %>% filter(GFP_RFP_b5<=q2 & GFP_RFP_b5>q1 & punc.RFP)
    data.q3= dtmp %>% filter(GFP_RFP_b5<=q3 & GFP_RFP_b5>q2 & punc.RFP)
    data.q4= dtmp %>% filter(GFP_RFP_b5<=q4 & GFP_RFP_b5>q3 & punc.RFP)
    data.q5= dtmp %>% filter(GFP_RFP_b5>q4 & punc.RFP)
    
    data=list(data.q1,data.q2,data.q3,data.q4,data.q5)
    
    frac.punc.q1 = nrow(data.q1)/sum(dtmp$GFP_RFP_b5<q1)
    frac.punc.q2 = nrow(data.q2)/sum(dtmp$GFP_RFP_b5<q2 & dtmp$GFP_RFP_b5>q1)
    frac.punc.q3 = nrow(data.q3)/sum(dtmp$GFP_RFP_b5<q3 & dtmp$GFP_RFP_b5>q2)
    frac.punc.q4 = nrow(data.q4)/sum(dtmp$GFP_RFP_b5<q4 & dtmp$GFP_RFP_b5>q3)
    frac.punc.q5 = nrow(data.q5)/sum(dtmp$GFP_RFP_b5>q4)
    #round the frac.punc list!!!#####
    frac.punc.list=list(round(frac.punc.q1,3),round(frac.punc.q2,3),round(frac.punc.q3,3),round(frac.punc.q4,3),round(frac.punc.q5,3))
    
    if (nrow(data.q1)>20 & nrow(data.q2)>20 & nrow(data.q3)>20 & nrow(data.q4)>20 & nrow(data.q5)>20){
      SUM=fig.data(n=5,FP="RFP")
      val.plot=val.plot.full(n=5,FP="RFP")
      RFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=5)
      csat=csat.plot(FP="RFP")
      RFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=5)
      frac.plot=frac.plot.full(FP="RFP")
      RFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20 & nrow(data.q2)>20 & nrow(data.q3)>20 & nrow(data.q4)>20){
      SUM=fig.data(n=4,FP="RFP")
      val.plot=val.plot.full(n=4,FP="RFP")
      RFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=4)
      csat=csat.plot(FP="RFP")
      RFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=4)
      frac.plot=frac.plot.full(FP="RFP")
      RFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20 & nrow(data.q2)>20 & nrow(data.q3)>20){
      SUM=fig.data(n=3,FP="RFP")
      val.plot=val.plot.full(n=3,FP="RFP")
      RFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=3)
      csat=csat.plot(FP="RFP")
      RFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=3)
      frac.plot=frac.plot.full(FP="RFP")
      RFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20 & nrow(data.q2)>20) {
      SUM=fig.data(n=2,FP="RFP")
      val.plot=val.plot.full(n=2,FP="RFP")
      RFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=2)
      csat=csat.plot(FP="RFP")
      RFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=2)
      frac.plot=frac.plot.full(FP="RFP")
      RFP.frac.plot.list[[I]]=frac.plot
      
    } else if (nrow(data.q1)>20) {
      SUM=fig.data(n=1,FP="RFP")
      val.plot=val.plot.full(n=1,FP="RFP")
      RFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=1)
      csat=csat.plot(FP="RFP")
      RFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=1)
      frac.plot=frac.plot.full(FP="RFP")
      RFP.frac.plot.list[[I]]=frac.plot
      
    } else {
      val.plot=val.plot.full(n=0,FP="RFP")
      RFP.valency.plot.list[[I]]=val.plot
      
      csat.data.df=csat.data(n=0)
      csat=csat.plot(FP="GFP")
      RFP.csat.plot.list[[I]]=csat
      
      frac.data.df=frac.data(n=0)
      frac.plot=frac.plot.full(FP="RFP")
      RFP.frac.plot.list[[I]]=frac.plot
    }
  }
  else{
    val.plot=val.plot.full(n=0,FP="RFP")
    RFP.valency.plot.list[[I]]=val.plot
    
    csat.data.df=csat.data(n=0)
    csat=csat.plot(FP="RFP")
    RFP.csat.plot.list[[I]]=csat
    
    frac.data.df=frac.data(n=0)
    frac.plot=frac.plot.full(FP="RFP")
    RFP.frac.plot.list[[I]]=frac.plot
    
  }
  #print(RFP.valency.plot.list[[I]])
  #print(RFP.csat.plot.list[[I]])
  #print(RFP.frac.plot.list[[I]])
}


panel.RFP=ggarrange(plotlist=RFP.valency.plot.list, ncol=1,nrow=22)
panel.sat.RFP=ggarrange(plotlist=RFP.csat.plot.list, ncol=1,nrow=22)
panel.frac.RFP=ggarrange(plotlist=RFP.frac.plot.list, ncol=1,nrow=22)

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/val_2_RFP_new.pdf",width=8,height=176)
print(panel.RFP)
dev.off()

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/csat_2_RFP_new.pdf",width=8,height=176)
print(panel.sat.RFP)
dev.off()

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt0/frac_2_RFP_new.pdf",width=8,height=176)
print(panel.frac.RFP)
dev.off()
