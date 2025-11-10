# Helper Functions: Analysis & Plotting

# Correlation

slope <- function(x, y){
  mean_x <- mean(x,na.rm=T)
  mean_y <- mean(y,na.rm=T)
  nom <- sum((x - mean_x)*(y-mean_y),na.rm=T)
  denom <- sum((x - mean_x)^2,na.rm=T)
  m <- nom / denom
  return(m)
}

# the slope formula is just
# covariance(x, y) / variance(x)
slope2 <- function(x, y){
  return(cov(x, y,use = 'pairwise', 'pearson')/var(x,na.rm = T))
}

intercept <- function(x, y, m){
  b <- mean(y,na.rm=T) - (m * mean(x,na.rm=T))
  return(b)
}
# Correlation with spearman method (ranks) and pairwise value
scor <- function(x,y,met='spearman',use='pairwise.complete.obs'){
  res = cor.test(x,y,method = met, use=use, exact=F) |>
    replace_name("r"="estimate","p"="p.value")
  res$pv = res$p.value
  res$p.value = ifelse(res$pv==0,"<1e-324" ,sprintf("%.1e",res$pv))
  return(res)
}
pearson <- function(X,Y){
  library(broom)
  res = cor.test(x = X, y=Y , method = "pearson",use='pairwise.complete',exact=F) |>
    broom::tidy() |>
    mutate( pv = p.value,
            p.value= ifelse(pv==0,"<1e-324" ,sprintf("%.1e",pv))) |>
    dplyr::rename(r=estimate,p=p.value)
  return(res)
}

spearman <- function(X,Y){
  library(broom)
  res = cor.test(x = X, y=Y , method = "spearman",use='pairwise.complete',exact=F) |>
    broom::tidy() |>
    mutate( pv = p.value,
            p.value= ifelse(pv==0,"<1e-324" ,sprintf("%.1e",pv))) |>
    dplyr::rename(r=estimate,p=p.value)
  return(res)
}

# Get spearman correlation parameters ready to plot
spearman.toplot = function(X,Y,rm.slope=T){
  s   = spearman(X,Y)
  s$slope = slope(X,Y)
  s$N = sum(complete.cases(X,Y))
  s$toshow = sprintf("%.3f\n%.3f\n%4s\n%4s",s$slope,s$r,s$p,s$N)
  if(rm.slope){
    s$toshow = sprintf("%.3f\n%4s\n%4s",s$r,s$p,s$N)
  }
  s$xmax = max(X,na.rm=T)
  s$ymax = max(Y,na.rm=T)
  s$xmin = min(X,na.rm=T)
  s$ymin = min(Y,na.rm=T)
  return(s)
}

pearson.toplot = function(X,Y,rm.slope=T){
  p   = pearson(X,Y)
  p$slope = slope(X,Y)
  p$N = sum(complete.cases(X,Y))
  p$toshow = sprintf("%.3f\n%.3f\n%4s\n%4s",p$slope,p$r,p$p,p$N)
  if(rm.slope){
    p$toshow = sprintf("%.3f\n%4s\n%4s",p$r,p$p,p$N)
  }
  p$xmax = max(X,na.rm=T)
  p$ymax = max(Y,na.rm=T)
  p$xmin = min(X,na.rm=T)
  p$ymin = min(Y,na.rm=T)
  return(p)
}

make_features_correlation= function(df_data){
  
  df_num = df_data %>%
    dplyr::select(where(~ is.numeric(.x))) %>%
    dplyr::select(-START,-END)
  
  cor_features = ggcorrplot::ggcorrplot(cor(df_num,use = 'pairwise'),
                                        outline.color = 'transparent',
                                        lab = T, lab_size = 1, lab_col='gray70',digits = 1,
                                        tl.cex = 8, tl.srt = 70)
  return(cor_features)
}

make_umap = function(df_data, K=10, seed=142, is_scaled=F){
  
  library(umap)
  library(ggalt)
  library(ggiraph)
  library(ggforce)
  set.seed(seed)
  umap.config = umap.defaults
  umap.config$n_neighbors = K
  umap.config$spread = 1
  umap.config$min_dist = 0.5
  umap.config$n_epochs = 500
  
  tag_neighbor = sprintf("-K%s",K)
  
  df_num = df_data %>%
    dplyr::select(where(~ is.numeric(.x))) %>%
    dplyr::select(-START,-END)
  
  df_info = df_data %>% dplyr::select( -colnames(df_num) )
  # Compute umap based on scaled features
  cat(sprintf("Compute umap with K=%s neighbors...\n",K))
  umap_data = umap::umap(d = as.matrix(df_num), seed = seed, config=umap.config)
  
  # Mark outliers on umap coordinates (sort distance and find elbow point)
  avgdist = umap_data$layout %>% dist %>% as.matrix() %>% colMeans()
  df_avgdist = tibble(x=rank(avgdist), y=avgdist)  %>%  mutate(idx = seq(avgdist)) %>% arrange(avgdist)
  elbow = pathviewr::find_curve_elbow(data_frame = df_avgdist[,1:2], export_type = 'all') 
  n_outliers = sum(avgdist > elbow['y'])
  cat(sprintf("%3d outliers average distance > %.2f\n",n_outliers,elbow['y']))
  
  pathviewr::find_curve_elbow(data_frame = df_avgdist[,1:2],plot_curve = T)
  elbow_plot = ggplot(df_avgdist) + 
    geom_point(aes(x=x,y=y)) +
    geom_vline(xintercept = elbow['x'],col='red') + 
    geom_text(x=elbow['x']-300,y=elbow['y']+1, col='red',
              label=sprintf('Elbow distance = %.2f\noutliers=%s',elbow['y'],n_outliers),
              hjust='inward',check_overlap = T) +
    ylab('UMAP average distance') + xlab('')
  
  df_umap_ = umap_data$layout %>% 
    magrittr::set_colnames(c('X1','X2')) %>%
    bind_cols(df_num) %>%
    bind_cols(df_info) %>%
    mutate(PROTEIN = ifelse(is.na(PROTEIN),AC,PROTEIN),
           pt_lab=paste0(PROTEIN," ",START,"-",END),
           outliers = avgdist > elbow['y'])
  
  x1.q = quantile(df_umap_$X1,seq(0,1,len=101))
  x2.q = quantile(df_umap_$X2,seq(0,1,len=101))
  print(summary(df_umap_[,c('X1','X2')]))
  
  df_umap = df_umap_ %>% filter(!outliers | from_atar)
  n_idr= n_distinct(df_info$IDR_id)
  n_idr_umap= n_distinct(df_umap$IDR_id)
  
  n_prot= n_distinct(df_info$AC)
  n_prot_umap= n_distinct(df_umap$AC)
  umap_criteria_text = sprintf("seed %s\n%s Neighbors\nscaled %s",seed,K,is_scaled)
  sample_size_text = sprintf("IDR=%s (%s)\nPROT=%s (%s)\noutliers=%s (D>%.2f)",n_idr_umap,n_idr,n_prot_umap,n_prot,n_outliers,elbow['y'])
  umap_coordinates_text = sprintf("X [%.1f,%.1f]\nY [%.1f,%.1f]",x1.q[1],last(x1.q),x2.q[1],last(x2.q))
  umap_atar = subset(df_umap,from_atar)
  # Plot the umap
  UMAP = ggplot(data=df_umap ,aes(x = X1,y = X2)) +
    # all idrs
    geom_point(size=1.5,shape=16,color='gray',alpha=1) +
    geom_text(data=NULL,label=sample_size_text, x=Inf, y=Inf, hjust='right',vjust='top',check_overlap = T, size=2) +
    geom_text(data=NULL,label=umap_criteria_text, x=-Inf, y=Inf, hjust='left',vjust='top',check_overlap = T, size=2) +
    geom_text(data=NULL,label=umap_coordinates_text, x=-Inf, y=-Inf, hjust='left',vjust='bottom',check_overlap = T, size=2) +
    
    # highlight atar idr
    geom_point(data=umap_atar, aes(color=PROTEIN), shape=16, size=3,alpha=0.9) +
    # annotate atar idr 
    ggrepel::geom_text_repel(data=umap_atar,aes(label=pt_lab,color=PROTEIN), size=3, fontface='italic',
                             max.overlaps=15, force=5, force_pull=0, seed=291222) +
    # graphical parameters
    labs(x = "UMAP1", y = "UMAP2", subtitle="")+
    see::scale_color_metro(palette = 'rainbow', discrete = T) +
    theme(legend.position="bottom") +
    ggeasy::easy_text_size(c("axis.title","axis.text.x", "axis.text.y"), size = 8)+
    ggeasy::easy_remove_legend() +
    ggeasy::easy_remove_x_axis('ticks') + ggeasy::easy_remove_y_axis('ticks') +
    theme(aspect.ratio = 1)
  
  # Make interactive points
  plot(UMAP)
  
  return(UMAP)
}
