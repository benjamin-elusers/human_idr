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
  "/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/BFP/"


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
colnames(platedef)=c("WELL","BF","RFP","GFP","BFP")
GFP.list.full=c("Ddx4","Dyrk3","FUS","GATA3","G3BP1","hnRNPA1",
                "HP1a",
                "HspB8","Nephrin","NPM1","RBM14","Syn1",
                "TAF15","TDP43","Valency","Syn2")
RFP.list.full=c("CycT1","Ddx3","FUS","GATA3","G3BP1","hnRNPA1",
                "HP1a","HspB8","Nephrin","NPM1","RBM14","Syn1",
                "TAF15","TDP43","Valency","Ddx4","DYRK3","Era","NUP98",
                "NUP98_trunc","UBQ2","Syn2")

GFPmono.list.full=c("CycT1_mono","FUS_mono","GATA3_mono","G3BP1_mono",
                    "hnRNPA1_mono","HP1a_mono","HspB8_mono",
                    "Nephrin_mono","NPM1_mono","RBM14_mono",
                    "TAF15_mono","TDP43_mono","Ddx4_mono","Era_mono",
                    "NUP98_mono","UBQ2_mono")

mut.list.full=c("hnRNPA1_WT","hnRNPA1_Aro+","hnRNPA1_Aro-",
                "hnRNPA1_Aro--","hnRNPA1_Aro_Perfect","hnRNPA1_Aro_Patchy",
                "TDP43_WT","TDP43_5StoD","TDP43_12StoD",
                "TDP43_A321V","TDP43_A326P","TDP43_320_330",
                "TDP43_331_343","TDP43_320_343")

###################################################################################

y.max <- c()
y.min <- c()
LLPS.reg.list <- list()
sat.qual <- c()
sat.val <- c()
frac.punc <- c()
GFP.plot.list <- list()
RFP.plot.list <- list()
GFPmono.plot.list <- list()
RFPmut.plot.list <- list()
GFPmut.plot.list <- list()



##### phase diagrams for RFP-tagged 2CLB-IDRs

for (I in c(1:length(RFP.list.full))){
  R = RFP.list.full[I]
  W = platedef$WELL[which(platedef$RFP==R)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("RFP-",R)
  punc = (dtmp$inRFPnfoci>0 & dtmp$f1_inRFP_toRFPmean>150)
  counts = sum(punc, na.rm=TRUE)
  frac.punc[I]=counts/ncell
  colors=c("salmon1","firebrick4")
  
  if (ncell>0 & counts>50){
    y.max[I] = quantile(dtmp$RFP_int_b5[punc],prob=0.90)
    y.min[I] =  quantile(dtmp$RFP_int_b5[punc],prob=0.20)
    punc.GOOD = punc & dtmp$RFP_int_b5<y.max[I] & dtmp$RFP_int_b5>y.min[I]
    LLPS.reg = lm(log(dtmp$RFP_int_b5[punc.GOOD])~log(dtmp$RFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] = LLPS.reg
    sat.val[I]=mean(dtmp$RFP_int_b5[punc.GOOD])
    sat.qual[I]=LLPS.reg$coefficients[2]
    plot.RFP=medmax.plot(FP="RFP",ST=1)
    RFP.plot.list[[I]]<-plot.RFP
  } else if (ncell>0 & counts>0) {
    plot.RFP=medmax.plot(FP="RFP",ST=2)
    RFP.plot.list[[I]]<-plot.RFP
  } else {
    plot.RFP=medmax.plot(FP="RFP",ST=3)
    RFP.plot.list[[I]]<-plot.RFP
  }
  print(RFP.plot.list[[I]])
}

panel.RFP=ggarrange(plotlist=RFP.plot.list, ncol=1,nrow=22) 
#  ggexport(filename = "/Users/atarg/Documents/weizmann/levy_lab/projects/llps_idr/imaging/2023_02_06_IDRt3/phase_diagram_GFPnew.pdf")

pdf.options(reset = TRUE, onefile = FALSE)
pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_16_hapIDRfull/phase_diagram_RFPsample_new.pdf",width=4,height=88)
print(panel.RFP)
dev.off()


#### phase diagrams for GFP-tagged 3IQ1-IDRs

for (I in c(1:length(GFP.list.full))){
  G = GFP.list.full[I]
  W = platedef$WELL[which(platedef$GFP==G)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("YFP-",G)
  punc = (dtmp$inGFPnfoci>0 & dtmp$f1_inGFP_toGFPmean>150)
  counts = sum(punc, na.rm=TRUE)
  frac.punc[I]=counts/ncell
  colors=c("lightgreen","darkgreen")
  
  if (ncell>0 & counts>50){
    y.max[I] = quantile(dtmp$GFP_int_b5[punc],prob=0.90)
    y.min[I] =  quantile(dtmp$GFP_int_b5[punc],prob=0.20)
    punc.GOOD = punc & dtmp$GFP_int_b5<y.max[I] & dtmp$GFP_int_b5>y.min[I]
    LLPS.reg = lm(log(dtmp$GFP_int_b5[punc.GOOD])~log(dtmp$GFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] = LLPS.reg
    sat.val[I]=mean(dtmp$GFP_int_b5[punc.GOOD])
    sat.qual[I]=LLPS.reg$coefficients[2]
    plot.GFP=medmax.plot(FP="GFP",ST=1)
    GFP.plot.list[[I]]<-plot.GFP
  } else if (ncell>0 & counts>0) {
    plot.GFP=medmax.plot(FP="GFP",ST=2)
    GFP.plot.list[[I]]<-plot.GFP
  } else {
    plot.GFP=medmax.plot(FP="GFP",ST=3)
    GFP.plot.list[[I]]<-plot.GFP
  }
  print(GFP.plot.list[[I]])
}

panel.GFP=ggarrange(plotlist=GFP.plot.list, ncol=1,nrow=16) 
#  ggexport(filename = "/Users/atarg/Documents/weizmann/levy_lab/projects/llps_idr/imaging/2023_02_06_IDRt3/phase_diagram_GFPnew.pdf")

pdf.options(reset = TRUE, onefile = FALSE)
pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_16_hapIDRfull/phase_diagram_GFPsample_new.pdf",width=4,height=64)
print(panel.GFP)
dev.off()


#### phase diagram for GFP-tagged monomeric IDRs

for (I in c(1:length(GFPmono.list.full))){
  G = GFPmono.list.full[I]
  W = platedef$WELL[which(platedef$GFP==G)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("YFP-",G)
  punc = (dtmp$inGFPnfoci>0 & dtmp$f1_inGFP_toGFPmean>150)
  counts = sum(punc, na.rm=TRUE)
  frac.punc[I]=counts/ncell
  colors=c("lightgreen","darkgreen")
  
  if (ncell>0 & counts>50){
    y.max[I] = quantile(dtmp$GFP_int_b5[punc],prob=0.90)
    y.min[I] =  quantile(dtmp$GFP_int_b5[punc],prob=0.20)
    punc.GOOD = punc & dtmp$GFP_int_b5<y.max[I] & dtmp$GFP_int_b5>y.min[I]
    LLPS.reg = lm(log(dtmp$GFP_int_b5[punc.GOOD])~log(dtmp$GFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] = LLPS.reg
    sat.val[I]=mean(dtmp$GFP_int_b5[punc.GOOD])
    sat.qual[I]=LLPS.reg$coefficients[2]
    plot.GFP=medmax.plot(FP="GFP",ST=1)
    GFPmono.plot.list[[I]]<-plot.GFP
  } else if (ncell>0 & counts>0) {
    plot.GFP=medmax.plot(FP="GFP",ST=2)
    GFPmono.plot.list[[I]]<-plot.GFP
  } else {
    plot.GFP=medmax.plot(FP="GFP",ST=3)
    GFPmono.plot.list[[I]]<-plot.GFP
  }
  print(GFPmono.plot.list[[I]])
}

panel.GFPmono=ggarrange(plotlist=GFPmono.plot.list, ncol=1,nrow=16) 
#  ggexport(filename = "/Users/atarg/Documents/weizmann/levy_lab/projects/llps_idr/imaging/2023_02_06_IDRt3/phase_diagram_GFPsample.pdf")

pdf.options(reset = TRUE, onefile = FALSE)
pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_16_hapIDRfull/phase_diagram_GFPmonosample_new.pdf",width=4,height=64)
print(panel.GFPmono)
dev.off()


##### phase diagrams for mutant GFP-tagged hnRNPA1 and TDP43:

for (I in c(1:length(mut.list.full))){
  G = mut.list.full[I]
  W = platedef$WELL[which(platedef$GFP==G)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("YFP-",G)
  punc = (dtmp$inGFPnfoci>0 & dtmp$f1_inGFP_toGFPmean>150)
  counts = sum(punc, na.rm=TRUE)
  frac.punc[I]=counts/ncell
  colors=c("lightgreen","darkgreen")
  
  if (ncell>0 & counts>50){
    y.max[I] = quantile(dtmp$GFP_int_b5[punc],prob=0.90)
    y.min[I] =  quantile(dtmp$GFP_int_b5[punc],prob=0.20)
    punc.GOOD = punc & dtmp$GFP_int_b5<y.max[I] & dtmp$GFP_int_b5>y.min[I]
    LLPS.reg = lm(log(dtmp$GFP_int_b5[punc.GOOD])~log(dtmp$GFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] = LLPS.reg
    sat.val[I]=mean(dtmp$GFP_int_b5[punc.GOOD])
    sat.qual[I]=LLPS.reg$coefficients[2]
    plot.GFP=medmax.plot(FP="GFP",ST=1)
    GFPmut.plot.list[[I]]<-plot.GFP
  } else if (ncell>0 & counts>0) {
    plot.GFP=medmax.plot(FP="GFP",ST=2)
    GFPmut.plot.list[[I]]<-plot.GFP
  } else {
    plot.GFP=medmax.plot(FP="GFP",ST=3)
    GFPmut.plot.list[[I]]<-plot.GFP
  }
  print(GFPmut.plot.list[[I]])
}

panel.GFPmut=ggarrange(plotlist=GFPmut.plot.list, ncol=1,nrow=14) 
#  ggexport(filename = "/Users/atarg/Documents/weizmann/levy_lab/projects/llps_idr/imaging/2023_02_06_IDRt3/phase_diagram_GFPnew.pdf")

pdf.options(reset = TRUE, onefile = FALSE)
pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_16_hapIDRfull/phase_diagram_GFPmutsample_new.pdf",width=4,height=56)
print(panel.GFPmut)
dev.off()


##### phase diagrams for mutant RFP-tagged hnRNPA1 and TDP43:

for (I in c(1:length(mut.list.full))){
  R = mut.list.full[I]
  W = platedef$WELL[which(platedef$RFP==R)]
  dtmp = data.KPlib[[W]]
  ncell=nrow(dtmp)
  strain = paste("RFP-",R)
  punc = (dtmp$inRFPnfoci>0 & dtmp$f1_inRFP_toRFPmean>150)
  counts = sum(punc, na.rm=TRUE)
  frac.punc[I]=counts/ncell
  colors=c("salmon1","firebrick4")
  
  if (ncell>0 & counts>50){
    y.max[I] = quantile(dtmp$RFP_int_b5[punc],prob=0.90)
    y.min[I] =  quantile(dtmp$RFP_int_b5[punc],prob=0.20)
    punc.GOOD = punc & dtmp$RFP_int_b5<y.max[I] & dtmp$RFP_int_b5>y.min[I]
    LLPS.reg = lm(log(dtmp$RFP_int_b5[punc.GOOD])~log(dtmp$RFP_int_b0[punc.GOOD]))
    LLPS.reg.list[[I]] = LLPS.reg
    sat.val[I]=mean(dtmp$RFP_int_b5[punc.GOOD])
    sat.qual[I]=LLPS.reg$coefficients[2]
    plot.RFP=medmax.plot(FP="RFP",ST=1)
    RFPmut.plot.list[[I]]<-plot.RFP
  } else if (ncell>0 & counts>0) {
    plot.RFP=medmax.plot(FP="RFP",ST=2)
    RFPmut.plot.list[[I]]<-plot.RFP
  } else {
    plot.RFP=medmax.plot(FP="RFP",ST=3)
    RFPmut.plot.list[[I]]<-plot.RFP
  }
  print(RFPmut.plot.list[[I]])
}

panel.RFPmut=ggarrange(plotlist=RFPmut.plot.list, ncol=1,nrow=14) 
#  ggexport(filename = "/Users/atarg/Documents/weizmann/levy_lab/projects/llps_idr/imaging/2023_02_06_IDRt3/phase_diagram_GFPnew.pdf")

pdf.options(reset = TRUE, onefile = FALSE)
pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_16_hapIDRfull/phase_diagram_RFPmutsample_new.pdf",width=4,height=56)
print(panel.RFPmut)
dev.off()