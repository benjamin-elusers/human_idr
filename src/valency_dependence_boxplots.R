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
  "/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_06_mutdip/"


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
GFP.list.full=c("FUS","hnRNPA1",
                "HspB8","RBM14",
                "TAF15","TDP43")
RFP.list.full=c("FUS","hnRNPA1",
                "HspB8","RBM14",
                "TAF15","TDP43","Ddx4",
                "DYRK3","Era","UBQ2")
GFPmut.list.full=c("hnRNPA1_WT","TDP43_WT","TDP43_12StoD",
                   "TDP43_5StoD","hnRNPA1_Aro+",
                   "hnRNPA1_Aro-","hnRNPA1_Aro_Perfect",
                   "hnRNPA1_Aro--","hnRNPA1_Aro_Patchy")



############################################################################

GFP.csat.plot.list <- list()
sat.val.GFP=c()
sat.qual.GFP=c()

for (I in seq_along(GFPmut.list.full)) {
  G <- GFPmut.list.full[I]
  R <- "Valency"
  W <- platedef$WELL[which(platedef$RFP == R & platedef$GFP == G)]
  dtmp <- data.KPlib[[W]]
  ncell <- nrow(dtmp)
  strain <- paste("GFP-", G)
  dtmp$GFP_RFP_b5 <- as.numeric(dtmp$RFP_int_mean)
  punc.GFP <- (dtmp$inGFPnfoci > 0 & dtmp$f1_inGFP_toGFPmean > 150)
  
  if (ncell > 0) {
    # Calculate quintiles for GFP_RFP_b5
    dtmp$Q <- as.factor(cut(dtmp$GFP_RFP_b5, breaks = quantile(dtmp$GFP_RFP_b5, probs = seq(0, 1, 0.20), na.rm = TRUE), labels = c(1:5)))
    
    # Filter dtmp based on punc.GFP
    dtmp <- dtmp[punc.GFP, ]
    
    # Ensure GFP_int_b5 is numeric
    dtmp$GFP_int_b5 <- as.numeric(dtmp$GFP_int_b5)
    
    # Calculate median RFP_int_b5 and count for each quintile
    summary_stats <- dtmp %>%
      group_by(Q) %>%
      summarise(median = round(median(GFP_int_b5),digits=1), count = n())
    
    # Generate ggplot boxplot
    plot <- ggplot(data = subset(dtmp, !is.na(Q)), aes(x = Q, y = GFP_int_b5, fill = Q)) +
      geom_boxplot() +
      geom_text(data = summary_stats, aes(label = paste("Median:", median, "\nCount:", count), y = 4500), vjust = -0.2) +
      geom_jitter(aes(color = Q), size = 2, width = 0.1) +
      scale_fill_manual(values = c("lightgoldenrod1", "orange1", "firebrick1", "darkred", "darkmagenta")) +
      scale_color_manual(values = c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")) +
      scale_y_continuous(limits = c(0, 5000)) +
      labs(title = paste("Boxplot for", G, "strain"),
           x = "Quintile",
           y = "GFP Intensity") +
      theme_minimal()
    print(plot)
    # Store the plot into the list
    GFP.csat.plot.list[[I]] <- plot
  }
}

panel.sat=ggarrange(plotlist=GFP.csat.plot.list, ncol=1,nrow=10)




pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2023_03_06_mutdip/csat_test_GFP.pdf",width=8,height=80)
print(panel.sat)
dev.off()


######################################################


RFP.csat.plot.list <- list()
sat.val.GFP=c()
sat.qual.GFP=c()
RFP.list.full=c("NPM1","Nephrin","HP1a","NUP98")

for (I in seq_along(RFP.list.full)) {
  G <- RFP.list.full[I]
  R <- "Valency"
  W <- platedef$WELL[which(platedef$RFP == G & platedef$GFP == R)]
  dtmp <- data.KPlib[[W]]
  ncell <- nrow(dtmp)
  strain <- paste("RFP-", G)
  dtmp$GFP_RFP_b5 <- as.numeric(dtmp$GFP_int_mean)
  punc.RFP <- (dtmp$inRFPnfoci > 0 & dtmp$f1_inRFP_toRFPmean > 150)
  
  if (ncell > 0) {
    # Calculate quintiles for GFP_RFP_b5
    dtmp$Q <- as.factor(cut(dtmp$GFP_RFP_b5, breaks = quantile(dtmp$GFP_RFP_b5, probs = seq(0, 1, 0.20), na.rm = TRUE), labels = c(1:5)))
    
    # Filter dtmp based on punc.GFP
    dtmp <- dtmp[punc.RFP, ]
    
    # Ensure RFP_int_b5 is numeric
    dtmp$RFP_int_b5 <- as.numeric(dtmp$RFP_int_b5)
    
    # Calculate median RFP_int_b5 and count for each quintile
    summary_stats <- dtmp %>%
      group_by(Q) %>%
      summarise(median = round(median(RFP_int_b5),digits=1), count = n())
    
    # Generate ggplot boxplot
    plot <- ggplot(data = subset(dtmp, !is.na(Q)), aes(x = Q, y = RFP_int_b5, fill = Q)) +
      geom_boxplot() +
      geom_text(data = summary_stats, aes(label = paste("Median:", median, "\nCount:", count), y = 4500), vjust = -0.2) +
      geom_jitter(aes(color = Q), size = 2, width = 0.1) +
      scale_fill_manual(values = c("lightgoldenrod1", "darkolivegreen2", "springgreen1", "springgreen3", "darkgreen")) +
      scale_color_manual(values = c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")) +
      scale_y_continuous(limits = c(0, 5000)) +
      labs(title = paste("Boxplot for", G, "strain"),
           x = "Quintile",
           y = "RFP Intensity") +
      theme_minimal()
    print(plot)
    # Store the plot into the list
    RFP.csat.plot.list[[I]] <- plot
  }
}

panel.sat=ggarrange(plotlist=RFP.csat.plot.list, ncol=1,nrow=10)




pdf("/Users/atarg/Documents/Levy_lab/IDR_LLPS/imaging/2022_09_14_IDRt3/csat_test_RFP.pdf",width=8,height=80)
print(panel.sat)
dev.off()




##############################################

