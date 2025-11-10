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
  "/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2024_06_11/"


#Read .RDS raw data results files:
output.dir=
  paste0(screen.dir, "results/")

data.KPlib=
  readRDS(paste0(output.dir,"raw_list.RDS"))[[1]]
design.KPlib=
  readRDS(paste0(output.dir,"design.RDS"))

#function to check which pictures (s_###) belong to a well in screen:
get.sNUM(design.KPlib,plate=1,well="B06")
dtmp=data.KPlib[["A06"]]
ncell=nrow(dtmp)
punc.GFP= (dtmp$inGFPnfoci>0 & dtmp$f1_inGFP_toGFPmean>150)
dtmp.punc=dtmp %>% filter(punc.GFP==TRUE)
dtmp.no.punc=dtmp %>% filter (punc.GFP==FALSE)
counts=sum(punc.GFP, na.rm = TRUE)
frac.punc=counts/ncell


ncell=nrow(dtmp)
punc.RFP= (dtmp$inRFPnfoci>0 & dtmp$f1_inRFP_toRFPmean>150)
dtmp.punc=dtmp %>% filter(punc.RFP==TRUE)
dtmp.no.punc=dtmp %>% filter (punc.RFP==FALSE)
counts=sum(punc.RFP, na.rm = TRUE)
frac.punc=counts/ncell

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

GFP.list.full=c("FUS","GATA3","G3BP1","hnRNPA1",
                "HP1a",
                "HspB8","Nephrin","NPM1","RBM14",
                #"Syn1",
                "TAF15","TDP43",
                "Ddx4","Era","NUP98","UBQ2")

RFP.list.full=c("FUS","GATA3","G3BP1","hnRNPA1",
                "HP1a","HspB8","Nephrin","NPM1","RBM14","Syn1",
                "TAF15","TDP43",
                "Ddx4","DYRK3","Era","NUP98",
                "UBQ2",
                "Valency")

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

GFP.ctrl.list.full=c("EDC3_35","PAB1_35","EDC3_50","PAB1_50","EDC3_35_starv","PAB1_35_starv","EDC3_50_starv","PAB1_50_starv")
RFP.ctrl.list.full=c("Ddx4","DYRK3","Era","FUS","hnRNPA1","HspB8","RBM14","TAF15","TDP43","UBQ2")

###############################################################################

# Distance score matrix for select strains:

# define plate:

D_score = matrix(NA, ncol=length(RFP.ctrl.list.full), nrow=length(GFP.ctrl.list.full))
# cellm=matrix(NA, ncol=length(RFP.list), nrow=length(GFP.list))
colnames(D_score)=RFP.ctrl.list.full
rownames(D_score)=GFP.ctrl.list.full

data.KPlib.mod=list()

for(I in c(1:length(GFP.ctrl.list.full))){
  
  for (J in c(1:length(RFP.ctrl.list.full))){
    
    G = GFP.ctrl.list.full[I]
    R = RFP.ctrl.list.full[J]
    W = platedef$WELL[which(platedef$GFP==G & platedef$RFP == R)]
    DAT = data.KPlib[[W]]
    ## compute stuff
    if (nrow(data.KPlib[[W]])>0){
      i=W
      data.KPlib.mod[[i]]=filter(data.KPlib[[i]], data.KPlib[[i]]$inGFPnfoci>=1 & data.KPlib[[i]]$inRFPnfoci>=1 & data.KPlib[[i]]$f1_inGFP_toGFPmean>150)
      if (nrow(data.KPlib.mod[[i]])>60){
        ##change cols to numeric:  
        data.KPlib.mod[[i]]$f1_inGFPx=
          as.numeric(data.KPlib.mod[[i]]$f1_inGFPx)
        data.KPlib.mod[[i]]$f1_inGFPy=
          as.numeric(data.KPlib.mod[[i]]$f1_inGFPy)
        data.KPlib.mod[[i]]$f1_inRFPx=
          as.numeric(data.KPlib.mod[[i]]$f1_inRFPx)
        data.KPlib.mod[[i]]$f1_inRFPy=
          as.numeric(data.KPlib.mod[[i]]$f1_inRFPy)
        data.KPlib.mod[[i]]$f2_inGFPx=
          as.numeric(data.KPlib.mod[[i]]$f2_inGFPx)
        data.KPlib.mod[[i]]$f2_inGFPy=
          as.numeric(data.KPlib.mod[[i]]$f2_inGFPy)
        data.KPlib.mod[[i]]$f2_inRFPx=
          as.numeric(data.KPlib.mod[[i]]$f2_inRFPx)
        data.KPlib.mod[[i]]$f2_inRFPy=
          as.numeric(data.KPlib.mod[[i]]$f2_inRFPy)
        
        data.KPlib.mod[[i]][data.KPlib.mod[[i]]$inGFPnfoci==1 & 
                              data.KPlib.mod[[i]]$inRFPnfoci==1,"f2_inRFPx"] = 100000
        data.KPlib.mod[[i]][data.KPlib.mod[[i]]$inGFPnfoci==1 
                            & data.KPlib.mod[[i]]$inRFPnfoci==1,"f2_inRFPy"] = 100000
        
        ## Calculate possible differences between x,y points:  
        data.KPlib.mod[[i]]$x_diff1=
          (data.KPlib.mod[[i]]$f1_inGFPx - data.KPlib.mod[[i]]$f1_inRFPx)^2
        data.KPlib.mod[[i]]$y_diff1=
          (data.KPlib.mod[[i]]$f1_inGFPy - data.KPlib.mod[[i]]$f1_inRFPy)^2
        data.KPlib.mod[[i]]$x_diff2=
          (data.KPlib.mod[[i]]$f2_inGFPx - data.KPlib.mod[[i]]$f2_inRFPx)^2
        data.KPlib.mod[[i]]$y_diff2=
          (data.KPlib.mod[[i]]$f2_inGFPy - data.KPlib.mod[[i]]$f2_inRFPy)^2
        data.KPlib.mod[[i]]$x_diff3=
          (data.KPlib.mod[[i]]$f1_inGFPx - data.KPlib.mod[[i]]$f2_inRFPx)^2
        data.KPlib.mod[[i]]$y_diff3=
          (data.KPlib.mod[[i]]$f1_inGFPy - data.KPlib.mod[[i]]$f2_inRFPy)^2
        data.KPlib.mod[[i]]$x_diff4=
          (data.KPlib.mod[[i]]$f2_inGFPx - data.KPlib.mod[[i]]$f1_inRFPx)^2
        data.KPlib.mod[[i]]$y_diff4=
          (data.KPlib.mod[[i]]$f2_inGFPy - data.KPlib.mod[[i]]$f1_inRFPy)^2
        
        ## Calculate possible distances between two points  
        data.KPlib.mod[[i]]$distance1=
          sqrt(data.KPlib.mod[[i]]$x_diff1 + data.KPlib.mod[[i]]$y_diff1)
        data.KPlib.mod[[i]]$distance2=
          sqrt(data.KPlib.mod[[i]]$x_diff2 + data.KPlib.mod[[i]]$y_diff2)
        data.KPlib.mod[[i]]$distance3=
          sqrt(data.KPlib.mod[[i]]$x_diff3 + data.KPlib.mod[[i]]$y_diff3)
        data.KPlib.mod[[i]]$distance4=
          sqrt(data.KPlib.mod[[i]]$x_diff4 + data.KPlib.mod[[i]]$y_diff4)
        
        ## Find minimum of distance between two foci  
        data.KPlib.mod[[i]]$distance=
          apply(data.KPlib.mod[[i]][,c("distance1","distance2","distance3","distance4")],1,min)
        data.KPlib.mod[[i]]$distance=data.KPlib.mod[[i]]$distance*0.108
        ## calculate mean distance per strain  
        D_score[I,J]=
          median(data.KPlib.mod[[i]]$distance)
      } else {
        D_score[I,J]=
          NA
      }
    } else {
      D_score[I,J]=
        NA
    }
  }
}

mode(D_score)="numeric"

distance_control_GFP = t(apply(D_score, 1,
                               function(x){
                                 div = x/D_score[nrow(D_score),]
                                 return(div)
                               } ))
distance_control_RFP = apply(D_score, 2,
                             function(x){
                               div = x/D_score[,(ncol(D_score))]
                               return(div)
                             } )
mode(distance_control_GFP)="numeric"
mode(distance_control_RFP)="numeric"

list.matrices=list(D_score,distance_control_GFP,distance_control_RFP)
main.titles=list("Foci Co-localization\nMedian Distance","Foci Co-localization\nGFP Normalized","Foci Co-localization\nRFP Normalized")

pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2024_06_11/pheatmaps_new50.pdf",width=16,height=9)
for (i in c(1:length(list.matrices))){
  rownames(list.matrices[[i]]) <- GFP.ctrl.list.full
  colnames(list.matrices[[i]]) <- RFP.ctrl.list.full
  pheatmap(list.matrices[[i]], 
           cluster_rows=FALSE, 
           cluster_cols=FALSE,
           border_color = "black",
           na_col = "gray97",
           breaks = c(0,0.4,0.6,1.5),
           color=c("lightgoldenrod1","orange1","firebrick1"),
           main=main.titles[i],
           display_numbers =round(list.matrices[[i]],2),
           fontsize=9)
}
dev.off()


