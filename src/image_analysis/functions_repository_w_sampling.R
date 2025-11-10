## FUNCTIONS

#Fluorescence correction functions:
correct.function <- function(list,GFPcol,RFPcol){
  for (i in c(1:length(list))){
    if (nrow(list[[i]])>0){
      list[[i]]=list[[i]][which(list[[i]][,GFPcol]<60000),]
      list[[i]]=list[[i]][which(list[[i]][,RFPcol]<60000),]
      list[[i]][,GFPcol]= list[[i]][,GFPcol]-(0.02*list[[i]][,RFPcol])
      list[[i]][which(list[[i]][,GFPcol]<1), GFPcol]=1
    }
  }
  return(list)
}

filter.saturated = function(list,GFPcol){
  for (i in c(1:length(list))){
    if (nrow(list[[i]])>0){
      list[[i]]=list[[i]][which(list[[i]][,GFPcol]<60000),]
    }
  }
  return(list)
}

#########################

#medmax.plot FUNCTION: generates scatterplot phase diagram
# of med/max values. Takes arguments of FP (either "RFP" or "GFP")
# and ST (which plot to create based on how many cells with
# condensate are found, ST=1 cells>50, ST=2 cells>0 and <50,
# ST=3 cells=0)
medmax.plot <- function(FP,ST) {
  if (FP=="GFP") {
    if (ST==1){
      dtmp.PUNC = dtmp %>% filter (punc==TRUE)
      dtmp.noPUNC= dtmp %>% filter (punc==FALSE)
      plot=
        ggplot()+
        geom_point(data=slice_sample(dtmp.noPUNC,n=(1000-(round(frac.punc[I],3)*1000))),
                   mapping=
                     aes(x=GFP_int_b0,y=GFP_int_b5),
                   colour=colors[1],
                   size=1,
                   stroke=1)+
        geom_point(data=slice_sample(dtmp.PUNC,n=(round(frac.punc[I],3)*1000)),
                   mapping=
                     aes(x=GFP_int_b0,y=GFP_int_b5),
                   colour=colors[2],
                   size=1,
                   stroke=1)+
        scale_x_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        geom_abline(slope=sat.qual[I],
                    intercept=LLPS.reg.list[[I]]$coefficients[1],
                    color=colors[2],
                    size=1)+
        geom_abline(slope=1,
                    intercept=0,
                    color="grey")+
        theme_classic()+
        labs(title=paste0(strain, 
                          "\nn = ",ncell,
                          "\np = ",(round(frac.punc[I],3)*100),
                          "\ns = ",(round(sat.qual[I],2)),
                          "\nCsat = ",(round(sat.val[I]))),
             x="Maximal GFP Intensity",
             y="Median GFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    } else if (ST==2) {
      dtmp.PUNC = dtmp %>% filter (punc==TRUE)
      dtmp.noPUNC= dtmp %>% filter (punc==FALSE)
      plot=
        ggplot()+
        geom_point(data=slice_sample(dtmp.noPUNC,n=(1000-(round(frac.punc[I],3)*1000))),
                   mapping=
                     aes(x=GFP_int_b0,y=GFP_int_b5),
                   colour=colors[1],
                   size=1,
                   stroke=1)+
        geom_point(data=slice_sample(dtmp.PUNC,n=(round(frac.punc[I],3)*1000)),
                   mapping=
                     aes(x=GFP_int_b0,y=GFP_int_b5),
                   colour=colors[2],
                   size=1,
                   stroke=1)+
        scale_x_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        geom_abline(slope=1,
                    intercept=0,
                    color="grey")+
        theme_classic()+
        labs(title=paste0(strain, 
                          "\nn = ",ncell,
                          "\np = ",(round(frac.punc[I],3)*100),
                          "\ns = undefined",
                          "\nCsat = undefined"),
             x="Maximal GFP Intensity",
             y="Median GFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    } else if (ST==3) {
      plot=
        ggplot()+
        geom_point(dtmp,
                   mapping=
                     aes(x=GFP_int_b0,y=GFP_int_b5),
                   colour=colors[1],
                   size=1,
                   stroke=1)+
        scale_x_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        geom_abline(slope=1,
                    intercept=0,
                    color="grey")+
        theme_classic()+
        labs(title=paste0(strain, 
                          "\nn = ",ncell,
                          "\np = ",(round(frac.punc[I],3)*100),
                          "\ns = undefined",
                          "\nCsat = undefined"),
             x="Maximal GFP Intensity",
             y="Median GFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    }
    
  } else if (FP=="RFP") {
    if (ST==1) {
      dtmp.PUNC = dtmp %>% filter (punc==TRUE)
      dtmp.noPUNC= dtmp %>% filter (punc==FALSE)
      plot=
        ggplot()+
        geom_point(data=slice_sample(dtmp.noPUNC,n=(1000-(round(frac.punc[I],3)*1000))),
                   mapping=
                     aes(x=RFP_int_b0,y=RFP_int_b5),
                   colour=colors[1],
                   size=1,
                   stroke=1)+
        geom_point(data=slice_sample(dtmp.PUNC,n=(round(frac.punc[I],3)*1000)),
                   mapping=
                     aes(x=RFP_int_b0,y=RFP_int_b5),
                   colour=colors[2],
                   size=1,
                   stroke=1)+
        scale_x_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        geom_abline(slope=sat.qual[I],
                    intercept=LLPS.reg.list[[I]]$coefficients[1],
                    color=colors[2],
                    size=1)+
        geom_abline(slope=1,
                    intercept=0,
                    color="grey")+
        theme_classic()+
        labs(title=paste0(strain, 
                          "\nn = ",ncell,
                          "\np = ",(round(frac.punc[I],3)*100),
                          "\ns = ",(round(sat.qual[I],2)),
                          "\nCsat = ",(round(sat.val[I]))),
             x="Maximal RFP Intensity",
             y="Median RFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    } else if (ST==2) {
      dtmp.PUNC = dtmp %>% filter (punc==TRUE)
      dtmp.noPUNC= dtmp %>% filter (punc==FALSE)
      plot=
        ggplot()+
        geom_point(data=slice_sample(dtmp.noPUNC,n=(1000-(round(frac.punc[I],3)*1000))),
                   mapping=
                     aes(x=RFP_int_b0,y=RFP_int_b5),
                   colour=colors[1],
                   size=1,
                   stroke=1)+
        geom_point(data=slice_sample(dtmp.PUNC,n=(round(frac.punc[I],3)*1000)),
                   mapping=
                     aes(x=RFP_int_b0,y=RFP_int_b5),
                   colour=colors[2],
                   size=1,
                   stroke=1)+
        scale_x_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        geom_abline(slope=1,
                    intercept=0,
                    color="grey")+
        theme_classic()+
        labs(title=paste0(strain, 
                          "\nn = ",ncell,
                          "\np = ",(round(frac.punc[I],3)*100),
                          "\ns = undefined",
                          "\nCsat = undefined"),
             x="Maximal RFP Intensity",
             y="Median RFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    } else if (ST==3) {
      plot=
        ggplot()+
        geom_point(dtmp,
                   mapping=
                     aes(x=RFP_int_b0,y=RFP_int_b5),
                   colour=colors[1],
                   size=1,
                   stroke=1)+
        scale_x_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(1,65000),breaks=c(1,10,10^2,10^3,10^4,10^5))+
        geom_abline(slope=1,
                    intercept=0,
                    color="grey")+
        theme_classic()+
        labs(title=paste0(strain, 
                          "\nn = ",ncell,
                          "\np = ",(round(frac.punc[I],3)*100),
                          "\ns = undefined",
                          "\nCsat = undefined"),
             x="Maximal RFP Intensity",
             y="Median RFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    }
  }
  return(plot)
}


#FUNCTION: generates the data necessary for valency plots - geom_point + abline
# stores them into a SUM list of lists. Takes arguments n corresponding
# to number of quintiles with cells in them for which linear regression can
# be calculated 

fig.data <- function(n,FP) {
  y.max<-list()
  y.min<-list()
  punc.GOOD.list<-list()
  LLPS.reg.list<-list()
  sat.qual<-list(0,0,0,0,0)
  sat.val<-list(0,0,0,0,0)
  stddev<-list(0,0,0,0,0)
  popdev<-list(0,0,0,0,0)
  Q <- list("0.2","0.4","0.6","0.8","1.0")
  geom.list<-list()
  abline.list<-list()
  
  if (FP=="GFP") {
    
    colors<-c("lightgoldenrod1","orange1","firebrick1","darkred","darkmagenta")

    for (i in c(1:n)){
      punc=(data[[i]]$inGFPnfoci>0 & data[[i]]$f1_inGFP_toGFPmean>150)
      print(punc)
      y.max[[i]]=quantile(data[[i]]$GFP_int_b5[punc],prob=0.90)
      y.min[[i]]=quantile(data[[i]]$GFP_int_b5[punc],prob=0.20)
      punc.GOOD=punc & data[[i]]$GFP_int_b5<y.max[[i]] & data[[i]]$GFP_int_b5>y.min[[i]]
      punc.GOOD.list[[i]]=punc.GOOD
      LLPS.reg= lm(log(data[[i]]$GFP_int_b5[punc.GOOD])~log(data[[i]]$GFP_int_b0[punc.GOOD]))
      LLPS.reg.list[[i]]=LLPS.reg
      sat.qual[[i]]=LLPS.reg$coefficients[2]
      sat.val[[i]]=mean(data[[i]]$GFP_int_b5[punc.GOOD])
      stddev[[i]]=sd(data[[i]]$GFP_int_b5[punc.GOOD])
      geom=geom_point(data[[i]],
                      mapping=aes(x=GFP_int_b0,
                                  y=GFP_int_b5),
                      colour=colors[i],
                      size=1.1,
                      stroke=1.1)
      abline=geom_abline(slope=LLPS.reg$coefficients[2],
                         intercept=LLPS.reg$coefficients[1],
                         color=colors[i],
                         linewidth=1.1)
      geom.list[[i]]<-geom
      abline.list[[i]]<-abline
      popdev[[i]]=sqrt((frac.punc.list[[i]])*(1-(frac.punc.list[[i]]))/nrow(data[[i]]))
    }
    
    SUM=list(y.max,y.min,punc.GOOD.list,LLPS.reg.list,sat.qual,sat.val,stddev,geom.list,abline.list,Q,popdev)
    
  } else if (FP=="RFP"){
    
    colors<-c("lightgoldenrod1","darkolivegreen2","springgreen1","springgreen3","darkgreen")
    
    for (i in c(1:n)){
      punc=(data[[i]]$inRFPnfoci>0 & data[[i]]$f1_inRFP_toRFPmean>150)
      print(punc)
      y.max[[i]]=quantile(data[[i]]$RFP_int_b5[punc],prob=0.90)
      y.min[[i]]=quantile(data[[i]]$RFP_int_b5[punc],prob=0.20)
      punc.GOOD=punc & data[[i]]$RFP_int_b5<y.max[[i]] & data[[i]]$RFP_int_b5>y.min[[i]]
      punc.GOOD.list[[i]]=punc.GOOD
      LLPS.reg= lm(log(data[[i]]$RFP_int_b5[punc.GOOD])~log(data[[i]]$RFP_int_b0[punc.GOOD]))
      LLPS.reg.list[[i]]=LLPS.reg
      sat.qual[[i]]=LLPS.reg$coefficients[2]
      sat.val[[i]]=mean(data[[i]]$RFP_int_b5[punc.GOOD])
      stddev[[i]]=sd(data[[i]]$RFP_int_b5[punc.GOOD])
      geom=geom_point(data[[i]],
                      mapping=aes(x=RFP_int_b0,
                                  y=RFP_int_b5),
                      colour=colors[i],
                      size=1.1,
                      stroke=1.1)
      abline=geom_abline(slope=LLPS.reg$coefficients[2],
                         intercept=LLPS.reg$coefficients[1],
                         color=colors[i],
                         linewidth=1.1)
      geom.list[[i]]<-geom
      abline.list[[i]]<-abline
      popdev[[i]]=sqrt((frac.punc.list[[i]])*(1-(frac.punc.list[[i]]))/nrow(data[[i]]))
    }
    
    SUM=list(y.max,y.min,punc.GOOD.list,LLPS.reg.list,sat.qual,sat.val,stddev,geom.list,abline.list,Q,popdev)
  }
  
  return(SUM)
}

# FUNCTION: generates the quintile barplot of % of cells with condensate:

frac.data <- function(n) {
  if (n>0){
    Q=unlist(SUM[10])
    frac=unlist(frac.punc.list)
    stddev=unlist(SUM[11])
    frac.data=data.frame(Q,frac,stddev)
    return(frac.data)
  } else {
    Q=c("0.2","0.4","0.6","0.8","1.0")
    frac=c(0,0,0,0,0)
    stddev=c(0,0,0,0,0)
    frac.data=data.frame(Q,frac,stddev)
    return(frac.data)
  }
}

frac.plot.full <- function(FP) {
  if (FP=="GFP") {
    colors=c("lightgoldenrod1","orange1","firebrick1","darkred","darkmagenta")
    plot=ggplot(frac.data.df,aes(y=frac,x=Q))+
      geom_bar(position='dodge2',
               stat='identity',
               aes(fill=Q))+
      geom_text(aes(label = frac), vjust = -0.2)+
      labs(title=paste0(strain,": by Quintile"),
           x=" ",
           y="p (%)")+
      scale_y_continuous(limits=c(0,1))+
      scale_fill_manual(values=colors)+
      geom_errorbar(aes(x=Q,ymin=frac-stddev,ymax=frac+stddev), colour="black",alpha=0.9,linewidth=0.5,width=0.2)+
      theme_classic()+
      theme(plot.title=element_text(size=9,
                                    lineheight=1,
                                    hjust=0.5),
            axis.title=element_text(size=9),
            axis.text=element_text(size=9))
  } else if (FP=="RFP") {
    colors=c("lightgoldenrod1","darkolivegreen2","springgreen1","springgreen3","darkgreen")
    plot=ggplot(frac.data.df,aes(y=frac,x=Q))+
      geom_bar(position='dodge2',
               stat='identity',
               aes(fill=Q))+
      geom_text(aes(label = frac), vjust = -0.2)+
      labs(title=paste0(strain,": by Quintile"),
           x=" ",
           y="p (%)")+
      scale_y_continuous(limits=c(0,1))+
      scale_fill_manual(values=colors)+
      geom_errorbar(aes(x=Q,ymin=frac-stddev,ymax=frac+stddev), colour="black",alpha=0.9,linewidth=0.5,width=0.2)+
      theme_classic()+
      theme(plot.title=element_text(size=9,
                                    lineheight=1,
                                    hjust=0.5),
            axis.title=element_text(size=9),
            axis.text=element_text(size=9))
  }
  return(plot)
}


# FUNCTION: generates the quintile barplot with stdev according to the number of available datasets
csat.data <- function(n){
  if (n>0){
  Q=unlist(SUM[10])
  sat.val=unlist(SUM[6])
  stddev=unlist(SUM[7])
  csat.data=data.frame(Q,sat.val,stddev)
  return(csat.data)
  } else {
  Q=c("0.2","0.4","0.6","0.8","1.0")
  sat.val=c(0,0,0,0,0)
  stddev=c(0,0,0,0,0)
  csat.data=data.frame(Q,sat.val,stddev)
  return(csat.data)
  }
}

csat.plot <- function(FP){
  if (FP=="GFP"){  colors=c("lightgoldenrod1","orange1","firebrick1","darkred","darkmagenta")
  plot=ggplot(csat.data.df,aes(y=sat.val,x=Q))+
    geom_bar(position='dodge2',
             stat='identity',
             aes(fill=Q))+
    geom_text(aes(label = round(sat.val)), vjust = -0.2)+
    labs(title=paste0(strain,": by Quintile"),
         x=" ",
         y="Csat (A.U.)")+
    scale_y_continuous(limits=c(0,5000))+
    scale_fill_manual(values=colors)+
    geom_errorbar(aes(x=Q,ymin=sat.val-stddev,ymax=sat.val+stddev), colour="black",alpha=0.9,linewidth=0.5,width=0.2)+
    theme_classic()+
    theme(plot.title=element_text(size=9,
                                  lineheight=1,
                                  hjust=0.5),
          axis.title=element_text(size=9),
          axis.text=element_text(size=9))
  } else if (FP=="RFP") {
    colors=c("lightgoldenrod1","darkolivegreen2","springgreen1","springgreen3","darkgreen")
    plot=ggplot(csat.data.df,aes(y=sat.val,x=Q))+
      geom_bar(position='dodge2',
               stat='identity',
               aes(fill=Q))+
      geom_text(aes(label = round(sat.val)), vjust = -0.2)+
      labs(title=paste0(strain,": by Quintile"),
           x=" ",
           y="Csat (A.U.)")+
      scale_y_continuous(limits=c(0,5000))+
      scale_fill_manual(values=colors)+
      geom_errorbar(aes(x=Q,ymin=sat.val-stddev,ymax=sat.val+stddev), colour="black",alpha=0.9,linewidth=0.5,width=0.2)+
      theme_classic()+
      theme(plot.title=element_text(size=9,
                                    lineheight=1,
                                    hjust=0.5),
            axis.title=element_text(size=9),
            axis.text=element_text(size=9))
  }
  return(plot)
}

# FUNCTION: generates valency scatter plot according to the number of available datasets
val.plot.full <- function(n,FP){
  if (FP=="GFP"){
    if (n>0) {
      plot=ggplot()+
        SUM[8]+
        scale_x_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        SUM[9]+
        theme_classic()+
        labs(title=paste0(strain, "\nTotal Cell Number = ",ncell,
                          "\nCsat 0.2 = ", (round(SUM[[6]][[1]],3)),
                          "\nCsat 0.4 = ", (round(SUM[[6]][[2]],3)),
                          "\nCsat 0.6 = ", (round(SUM[[6]][[3]],3)),
                          "\nCsat 0.8 = ", (round(SUM[[6]][[4]],3)),
                          "\nCsat 1.0 = ", (round(SUM[[6]][[5]],3)),
                          "\n% 0.2 = ", (round(frac.punc.q1,3)*100),
                          "\n% 0.4 = ", (round(frac.punc.q2,3)*100),
                          "\n% 0.6 = ", (round(frac.punc.q3,3)*100),
                          "\n% 0.8 = ", (round(frac.punc.q4,3)*100),
                          "\n% 1.0 = ", (round(frac.punc.q5,3)*100)),
             x="Maximal GFP Intensity",
             y="Median GFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    } else if (n==0) {
      plot=ggplot()+
        geom_blank(mapping=NULL,
                   data=NULL)+
        scale_x_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        theme_classic()+
        labs(title=paste0(strain, "\nTotal Cell Number = ",ncell,
                          "\nCsat 0.2 = undefined",
                          "\nCsat 0.4 = undefined",
                          "\nCsat 0.6 = undefined",
                          "\nCsat 0.8 = undefined",
                          "\nCsat 1.0 = undefined",
                          "\n% 0.2 = undefined",
                          "\n% 0.4 = undefined",
                          "\n% 0.6 = undefined",
                          "\n% 0.8 = undefined",
                          "\n% 1.0 = undefined"),
             x="Maximal GFP Intensity",
             y="Median GFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    }
  } else if (FP=="RFP") {
    if (n>0) {
      plot=ggplot()+
        SUM[8]+
        scale_x_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        SUM[9]+
        theme_classic()+
        labs(title=paste0(strain, "\nTotal Cell Number = ",ncell,
                          "\nCsat 0.2 = ", (round(SUM[[6]][[1]],3)),
                          "\nCsat 0.4 = ", (round(SUM[[6]][[2]],3)),
                          "\nCsat 0.6 = ", (round(SUM[[6]][[3]],3)),
                          "\nCsat 0.8 = ", (round(SUM[[6]][[4]],3)),
                          "\nCsat 1.0 = ", (round(SUM[[6]][[5]],3)),
                          "\n% 0.2 = ", (round(frac.punc.q1,3)*100),
                          "\n% 0.4 = ", (round(frac.punc.q2,3)*100),
                          "\n% 0.6 = ", (round(frac.punc.q3,3)*100),
                          "\n% 0.8 = ", (round(frac.punc.q4,3)*100),
                          "\n% 1.0 = ", (round(frac.punc.q5,3)*100)),
             x="Maximal RFP Intensity",
             y="Median RFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    } else if (n==0) {
      plot=ggplot()+
        geom_blank(mapping=NULL,
                   data=NULL)+
        scale_x_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        scale_y_continuous(trans='log',limits =c(10,65000),breaks=c(10,10^2,10^3,10^4,10^5))+
        theme_classic()+
        labs(title=paste0(strain, "\nTotal Cell Number = ",ncell,
                          "\nCsat 0.2 = undefined",
                          "\nCsat 0.4 = undefined",
                          "\nCsat 0.6 = undefined",
                          "\nCsat 0.8 = undefined",
                          "\nCsat 1.0 = undefined",
                          "\n% 0.2 = undefined",
                          "\n% 0.4 = undefined",
                          "\n% 0.6 = undefined",
                          "\n% 0.8 = undefined",
                          "\n% 1.0 = undefined"),
             x="Maximal RFP Intensity",
             y="Median RFP Intensity")+
        theme(plot.title=element_text(size=9,
                                      lineheight=1,
                                      hjust=0.5),
              axis.title=element_text(size=9),
              axis.text=element_text(size=9))
    }
  }
  return(plot)
}

