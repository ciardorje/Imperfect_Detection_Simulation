rm(list=ls()); gc()

library(pacman)
p_load(ggplot2, ggpubr, ragg)

setwd()

#-----Plot Performance Across All Repetitions-----
load('Estimator_Performance_Results.RData')

#Plotting helpers
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
data_summary <- function(x) {
  
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
  
}

metrics <- c('Mean_Bias', 
             'Mean_Accuracy', 
             #'Mean_Coverage', 
             'Precision')

#####Plot function#####
plotty <- function(x, estimators, metrics, scenarios){
  
  #Plot holders
  sub_plot_list <- vector('list', length(metrics)) 
  names(sub_plot_list) <- metrics
  x_plots <- vector('list', length(scenarios))
  names(x_plots) <- scenarios
  x_plots <- lapply(x_plots, function(x) x_plots[x] <- sub_plot_list)
  
  #Plot labels
  nsamples_labels <- c(`2` = '2 Sampling Events', 
                       `3` = '3 Sampling Events',
                       `5` = '5 Sampling Events')
  
  for(s in scenarios){
    for(m in metrics){     
      
      #Isolate plot data
      plot_data <- x[x$Detection_Scenario == s, c('Estimator', 'Sampling_Events', m)]
      names(plot_data)[3] <- 'Metric'
      
      #Create base plot
      x_plots[[s]][[m]] <- 
        ggplot(plot_data,
               aes(y = Metric,
                   x = Estimator)) +
        geom_violin(alpha = 0.7, trim = F, size = 1, bw = 0.04, 
                    colour = '#DB0F14',  
                    fill = '#DB0F14') +
        stat_summary(fun.data = data_summary,
                     geom = "pointrange", color = "black", size = 0.3) +
        theme_bw() +
        ylab(sub('_', ' ', m)) +
        xlab(gsub('_', ' ', estimators)) +
        theme(legend.position = 'none',
              plot.margin = unit(c(3.5, 7.5, 3.5, 7.5), 'pt'),
              title = element_text(face = 'bold', size = 9),
              axis.title.x = element_blank(),
              axis.text.x = element_text(face = 'bold'), 
              axis.title.y = element_text(vjust = 2.25),
              axis.text.y = element_text(face = 'bold'))
      
      #Facet Plots and set sppropriate axes for each metric
      if(m == 'Mean_Bias'){
        
        x_plots[[s]][[m]] <-
          x_plots[[s]][[m]] +
          scale_y_continuous(limits = c(floor_dec(quantile(x$Mean_Bias, 0.025)),
                                        ceiling_dec(quantile(x$Mean_Bias, 0.975))),
                             breaks = round(seq(floor_dec(quantile(x$Mean_Bias, 0.025)),
                                                ceiling_dec(quantile(x$Mean_Bias, 0.975)),.4), 1),
                             labels = round(seq(floor_dec(quantile(x$Mean_Bias, 0.025)),
                                                ceiling_dec(quantile(x$Mean_Bias, 0.975)),.4), 1)) +
          geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
          facet_grid(~Sampling_Events, 
                     labeller = as_labeller(nsamples_labels)) +
          theme(strip.text = element_text(face = 'bold', size = 9.5)) 
        
      } else if(m == 'Mean_Accuracy'){
        
        x_plots[[s]][[m]] <-
          x_plots[[s]][[m]] +
          scale_y_continuous(limits = c(0, ceiling_dec(quantile(x$Mean_Accuracy, 0.975))),
                             breaks = seq(0,ceiling_dec(quantile(x$Mean_Accuracy, 0.975)),.4),
                             labels = seq(0,ceiling_dec(quantile(x$Mean_Accuracy, 0.975)),.4)) +
          facet_grid(~Sampling_Events) +
          theme(strip.text = element_blank()) 
        
        
      } else if(m == 'Mean_Coverage'){
        
        x_plots[[s]][[m]] <-
          x_plots[[s]][[m]] +
          scale_y_continuous(limits = c(0,1),
                             breaks = seq(0,1,.2),
                             labels = seq(0,1,.2)) +
          facet_grid(~Sampling_Events) +
          theme(strip.text = element_blank())
        
        
      } else if(m == 'Precision'){
        
        x_plots[[s]][[m]] <-
          x_plots[[s]][[m]] +
          scale_y_continuous(limits = c(0, ceiling_dec(quantile(x$Precision, 0.975))),
                             breaks = seq(0, ceiling_dec(quantile(x$Precision, 0.975)),.2),
                             labels = seq(0, ceiling_dec(quantile(x$Precision, 0.975)),.2)) +
          facet_grid(~Sampling_Events) +
          theme(strip.text = element_blank()) 
        
      }
    }
  }
  
  #Combine metric plots for each detection scenario
  combined_x_plots <- vector('list', length(scenarios))
  for(i in 1:6){
    
    combined_x_plots[[i]] <- ggarrange(x_plots[[i]][[1]], x_plots[[i]][[2]],
                                       x_plots[[i]][[3]], #x_plots[[i]][[4]],
                                       ncol = 1, nrow = length(metrics),
                                       heights = c(1.2,1,1#,1
                                       ))
    combined_x_plots[[i]]$theme$plot.margin <- unit(c(22.5, 3, 17, 3), 'pt')
    
  }
  
  return(combined_x_plots)
  
}


#####Site-Level Richness#####
siteR <- all_ests$Richness_Performance
siteR$Sampling_Events[siteR$Samples == 'Single'] <- 1  
siteR <- siteR[siteR$Sampling_Events > 1,] 
siteR <- siteR[,-c(6:14,20)] %>% distinct()
siteR <- siteR[siteR$Estimator == 'Observed' | siteR$Non_NA_Site >= 0.8,]
siteR$Precision[is.na(siteR$Precision)] <- 0
siteR <- siteR[!siteR$Estimator %in% c('ACE', 'Jack2'),]

#Factorise estimators to help plot discretely 
estimators <- c('Observed', 'iChao1', 'HMSOM', 'HMSOM_Predicted')
siteR$Estimator <- factor(siteR$Estimator, levels = estimators)

#Plot
combined_SiteR_plots <- plotty(siteR, 
                               estimators = estimators,
                               metrics = metrics,
                               scenarios = 1:6)

#Combine all plots and save as png
SiteR_full_plot <- fs::path('./Figures/',  "Site_Level_Richness.png")
agg_png(SiteR_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_SiteR_plots[[1]], combined_SiteR_plots[[4]],
          combined_SiteR_plots[[2]], combined_SiteR_plots[[5]],
          combined_SiteR_plots[[3]], combined_SiteR_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(SiteR_full_plot)

######Site-Level Beta#####
siteB <- all_ests$Beta_Performance
siteB$Sampling_Events[siteB$Samples == 'Single'] <- 1  
siteB <- siteB[siteB$Sampling_Events > 1,] 
siteB <- siteB[,-c(3:6,11:14,20)] %>% distinct()
siteB <- siteB[siteB$Estimator == 'Observed' | siteB$Non_NA_Pair >= 0.8,]
siteB$Precision[is.na(siteB$Precision)] <- 0

b_estimators <- c('Observed', 'iChao1', 'HMSOM', 'HMSOM_Predicted')
siteB$Estimator <- factor(siteB$Estimator, levels = b_estimators)

combined_SiteB_plots <- plotty(siteB, 
                               estimators = b_estimators,
                               metrics = metrics,
                               scenarios = 1:6)

SiteB_full_plot <- fs::path('./Figures/',  "Site_Level_Beta.png")
agg_png(SiteB_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_SiteB_plots[[1]], combined_SiteB_plots[[4]],
          combined_SiteB_plots[[2]], combined_SiteB_plots[[5]],
          combined_SiteB_plots[[3]], combined_SiteB_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1, vjust = 1,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(SiteB_full_plot)

#####SAR Slope#####
SARz <- all_ests$SAR_Performance
SARz <- separate(SARz, Estimator, into = c('Estimator', 'Samples'), sep = '\\.')
SARz <- SARz[is.na(SARz$Samples) | SARz$Samples != 'Single',] 
SARz <- SARz[SARz$Estimator == 'Observed' | SARz$Type == 'LM',]
SARz <- SARz[,c(1, 24:27, 34, 38:39)] %>% distinct()
names(SARz)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')
SARz <- SARz[SARz$Estimator == 'Observed' | SARz$Non_NA_Repeats >= 0.8,]
SARz <- SARz[!SARz$Estimator %in% c('ACE', 'Jack2'),]

SARz$Estimator <- factor(SARz$Estimator, levels = estimators)

combined_SARz_plots <- plotty(SARz, 
                              estimators = estimators,
                              metrics = metrics,
                              scenarios = 1:6)

SARz_full_plot <- fs::path('./Figures/',  "SAR_Slope.png")
agg_png(SARz_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_SARz_plots[[1]], combined_SARz_plots[[4]],
          combined_SARz_plots[[2]], combined_SARz_plots[[5]],
          combined_SARz_plots[[3]], combined_SARz_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(SARz_full_plot)


#####Beta Model Slope#####
Beta_z <- all_ests$BetaMod_Performance
Beta_z <- separate(Beta_z, Estimator, into = c('Estimator', 'Samples'), sep = '\\.')
Beta_z <- Beta_z[is.na(Beta_z$Samples) | Beta_z$Samples != 'Single',] 
Beta_z <- Beta_z[Beta_z$Estimator == 'Observed' | Beta_z$Type == 'LM',]
Beta_z$b_Precision[is.na(Beta_z$b_Precision)] <- 0
Beta_z <- Beta_z[,c(1, 24:27, 33, 35, 39)] %>% distinct()
names(Beta_z)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')
Beta_z <- Beta_z[Beta_z$Estimator == 'Observed' | Beta_z$Non_NA_Repeats >= 0.8,]
Beta_z <- Beta_z[!is.na(Beta_z$Mean_Bias),]

Beta_z$Estimator <- factor(Beta_z$Estimator, levels = b_estimators)

combined_Beta_z_plots <- plotty(Beta_z, 
                                estimators = b_estimators,
                                metrics = metrics,
                                scenarios = 1:6)

Betaz_full_plot <- fs::path('./Figures/',  "Beta_Model_Slope.png")
agg_png(Betaz_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_Beta_z_plots[[1]], combined_Beta_z_plots[[4]],
          combined_Beta_z_plots[[2]], combined_Beta_z_plots[[5]],
          combined_Beta_z_plots[[3]], combined_Beta_z_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(Betaz_full_plot)

#####SAR Intercept#####
SARc <- all_ests$SAR_Performance
SARc <- separate(SARc, Estimator, into = c('Estimator', 'Samples'), sep = '\\.')
SARc <- SARc[is.na(SARc$Samples) | SARc$Samples != 'Single',] 
SARc <- SARc[SARc$Estimator == 'Observed' | SARc$Type == 'LM',]
SARc <- SARc[,c(1, 28:31, 34, 38:39)] %>% distinct()
names(SARc)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')
SARc <- SARc[SARc$Estimator == 'Observed' | SARc$Non_NA_Repeats >= 0.8,]
SARc <- SARc[!SARc$Estimator %in% c('ACE', 'Jack2'),]

SARc$Estimator <- factor(SARc$Estimator, levels = estimators)

combined_SARc_plots <- plotty(SARc, 
                              estimators = estimators,
                              metrics = metrics,
                              scenarios = 1:6)

SARc_full_plot <- fs::path('./Figures/',  "SAR_Intercept.png")
agg_png(SARc_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_SARc_plots[[1]], combined_SARc_plots[[4]],
          combined_SARc_plots[[2]], combined_SARc_plots[[5]],
          combined_SARc_plots[[3]], combined_SARc_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(SARc_full_plot)

#####Beta Model Intercept#####
Betac <- all_ests$BetaMod_Performance
Betac <- separate(Betac, Estimator, into = c('Estimator', 'Samples'), sep = '\\.')
Betac <- Betac[is.na(Betac$Samples) | Betac$Samples != 'Single',] 
Betac <- Betac[Betac$Estimator == 'Observed' | Betac$Type == 'LM',]
Betac <- Betac[,c(1, 28:31, 35, 39, 33)] %>% distinct()
names(Betac)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')
Betac <- Betac[Betac$Estimator == 'Observed' | Betac$Non_NA_Repeats >= 0.8,]
Betac <- Betac[!is.na(Betac$Mean_Bias),]

Betac$Estimator <- factor(Betac$Estimator, levels = b_estimators)

combined_Beta_c_plots <- plotty(Betac, 
                                estimators = b_estimators,
                                metrics = metrics,
                                scenarios = 1:6)

Betac_full_plot <- fs::path('./Figures/',  "Beta_Intercept.png")
agg_png(Betac_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_Beta_c_plots[[1]], combined_Beta_c_plots[[4]],
          combined_Beta_c_plots[[2]], combined_Beta_c_plots[[5]],
          combined_Beta_c_plots[[3]], combined_Beta_c_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'),  
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(Betac_full_plot)


#####Combined intercept and slope accuracy (overall model fit)#####
#SARs
SARov <- all_ests$SAR_Performance
SARov <- separate(SARov, Estimator, into = c('Estimator', 'Samples'), sep = '\\.')
SARov <- SARov[is.na(SARov$Samples) | SARov$Samples != 'Single',] 
SARov <- SARov[SARov$Estimator == 'Observed' | SARov$Type == 'LM',]
SARov <- SARov[,c(1, 32, 34, 38:39)] %>% distinct()
SARov <- SARov[SARov$Estimator == 'Observed' | SARov$Non_NA_Repeats >= 0.8,]
SARov <- SARov[!SARov$Estimator %in% c('ACE', 'Jack2'),]
names(SARov)[2] <- 'Mean_Accuracy'

SARov$Estimator <- factor(SARov$Estimator, levels = estimators)

SARov_plots <- plotty(SARov,
                      estimators = estimators,
                      metrics = 'Mean_Accuracy',
                      scenarios = 1:6)

SARov_full_plot <- fs::path('./Figures/',  "SAR_Overall_Accuracy.png")
agg_png(SARov_full_plot, width = 600, height = 400, units = "mm", res = 500) 

ggarrange(SARov_plots[[1]], SARov_plots[[4]],
          SARov_plots[[2]], SARov_plots[[5]],
          SARov_plots[[3]], SARov_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1,
          ncol = 2, nrow = 3,
          heights = c(1.05, 1, 1))


invisible(dev.off())
knitr::include_graphics(SARov_full_plot)

#Beta models
Bov <- all_ests$BetaMod_Performance
Bov <- separate(Bov, Estimator, into = c('Estimator', 'Samples'), sep = '\\.')
Bov <- Bov[is.na(Bov$Samples) | Bov$Samples != 'Single',] 
Bov <- Bov[Bov$Estimator == 'Observed' | Bov$Type == 'LM',]
Bov <- Bov[,c(1, 32, 33, 35, 39)] %>% distinct()
Bov <- Bov[Bov$Estimator == 'Observed' | Bov$Non_NA_Repeats >= 0.8,]

Bov$Estimator <- factor(Bov$Estimator, levels = b_estimators)

Bov_plots <- plotty(Bov,
                    estimators = b_estimators,
                    metrics = 'Mean_Accuracy',
                    scenarios = 1:6)

Bov_full_plot <- fs::path('./Figures/',  "Beta_Overall_Accuracy.png")
agg_png(Bov_full_plot, width = 600, height = 400, units = "mm", res = 500) 

ggarrange(Bov_plots[[1]], Bov_plots[[4]],
          Bov_plots[[2]], Bov_plots[[5]],
          Bov_plots[[3]], Bov_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1,
          ncol = 2, nrow = 3,
          heights = c(1.05, 1, 1))


invisible(dev.off())
knitr::include_graphics(Bov_full_plot)

######Baselga plots######

#Site-level 
#Turnover
siteBaselga <- all_ests$Baselga_Performance

siteT <- siteBaselga[,c(1:3, 22:25, 30:31, 35)] %>% distinct()
siteN <- siteBaselga[,c(1:3,26:29,30:31,35)] %>% distinct()
names(siteT)[4:7] <- names(siteN)[4:7] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')

baselga_estimators <- c('Observed', 'HMSOM', 'HMSOM_Predicted')
siteT$Estimator <- factor(siteT$Estimator, levels = baselga_estimators)
siteN$Estimator <- factor(siteN$Estimator, levels = baselga_estimators)

siteT$Precision[is.na(siteT$Precision)] <- 0
siteN$Precision[is.na(siteN$Precision)] <- 0

combined_siteT_plots <- plotty(siteT,
                               estimators = baselga_estimators,
                               metrics = metrics,
                               scenarios = 1:6)

siteT_full_plot <- fs::path('./Figures/',  "Site_Level_Turnover.png")
agg_png(siteT_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_siteT_plots[[1]], combined_siteT_plots[[4]],
          combined_siteT_plots[[2]], combined_siteT_plots[[5]],
          combined_siteT_plots[[3]], combined_siteT_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1, vjust = 1,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(siteT_full_plot)

#Nestedness
combined_siteN_plots <- plotty(siteN,
                               estimators = baselga_estimators,
                               metrics = metrics,
                               scenarios = 1:6)

siteN_full_plot <- fs::path('./Figures/',  "Site_Level_Nestedness.png")
agg_png(siteN_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_siteN_plots[[1]], combined_siteN_plots[[4]],
          combined_siteN_plots[[2]], combined_siteN_plots[[5]],
          combined_siteN_plots[[3]], combined_siteN_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'),  
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1, vjust = 1,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(siteN_full_plot)

#Baselga models
#Slopes
Baselga_b <- all_ests$BaselgaMod_Performance
Baselga_b$b_Precision[is.na(Baselga_b$b_Precision)] <- 0

Turn_b <- Baselga_b[Baselga_b$Variable == 'Turnover', c(1:3, 21:24, 30, 31, 35)] %>% distinct()
Nest_b <- Baselga_b[Baselga_b$Variable == 'Nestedness', c(1:3, 25:28, 30, 31, 35)] %>% distinct()

Turn_b <- Turn_b[Turn_b$Estimator == 'Observed' | Turn_b$Type == 'LM',]
Nest_b <- Nest_b[Nest_b$Estimator == 'Observed' | Nest_b$Type == 'LM',]

names(Nest_b)[4:7] <- names(Turn_b)[4:7] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')

Turn_b$Estimator <- factor(Turn_b$Estimator, levels = baselga_estimators)
Nest_b$Estimator <- factor(Nest_b$Estimator, levels = baselga_estimators)

#Turnover
combined_Turn_b_plots <- plotty(Turn_b,
                                estimators = baselga_estimators,
                                metrics = metrics,
                                scenarios = 1:6)

Turnb_full_plot <- fs::path('./Figures/',  "Turnover_Model_Slope.png")
agg_png(Turnb_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_Turn_b_plots[[1]], combined_Turn_b_plots[[4]],
          combined_Turn_b_plots[[2]], combined_Turn_b_plots[[5]],
          combined_Turn_b_plots[[3]], combined_Turn_b_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(Turnb_full_plot)

#Nestedness
combined_Nest_b_plots <- plotty(Nest_b,
                                estimators = baselga_estimators,
                                metrics = metrics,
                                scenarios = 1:6)


Nest_b_full_plot <- fs::path('./Figures/',  "Nestedness_Model_Slope.png")
agg_png(Nest_b_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_Nest_b_plots[[1]], combined_Nest_b_plots[[4]],
          combined_Nest_b_plots[[2]], combined_Nest_b_plots[[5]],
          combined_Nest_b_plots[[3]], combined_Nest_b_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(Nest_b_full_plot)

#Intercepts
Baselga_a <- all_ests$BaselgaMod_Performance

Baselga_a$a_Precision[is.na(Baselga_a$a_Precision)] <- 0
Turn_a <- Baselga_a[Baselga_a$Variable == 'Turnover', c(1:3, 25:28, 30, 31, 35)] %>% distinct()
Nest_a <- Baselga_a[Baselga_a$Variable == 'Nestedness', c(1:3, 25:28, 30, 31, 35)] %>% distinct()

Turn_a <- Turn_a[Turn_a$Estimator == 'Observed' | Turn_a$Type == 'LM',]
Nest_a <- Nest_a[Nest_a$Estimator == 'Observed' | Nest_a$Type == 'LM',]

names(Turn_a)[4:7] <- names(Nest_a)[4:7] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')

Turn_a$Estimator <- factor(Turn_a$Estimator, levels = baselga_estimators)
Nest_a$Estimator <- factor(Nest_a$Estimator, levels = baselga_estimators)

combined_Turn_a_plots <- plotty(Turn_a,
                                estimators = baselga_estimators,
                                metrics = metrics,
                                scenarios = 1:6)

Turn_a_full_plot <- fs::path('./Figures/',  "Turnover_Intercept.png")
agg_png(Turn_a_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_Turn_a_plots[[1]], combined_Turn_a_plots[[4]],
          combined_Turn_a_plots[[2]], combined_Turn_a_plots[[5]],
          combined_Turn_a_plots[[3]], combined_Turn_a_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(Turn_a_full_plot)

#Nestedness
combined_Nest_a_plots <- plotty(Nest_a,
                                estimators = baselga_estimators,
                                metrics = metrics,
                                scenarios = 1:6)

Nest_a_full_plot <- fs::path('./Figures/',  "Nestedness_Intercept.png")
agg_png(Nest_a_full_plot, width = 600, height = 600, units = "mm", res = 500) 

ggarrange(combined_Nest_a_plots[[1]], combined_Nest_a_plots[[4]],
          combined_Nest_a_plots[[2]], combined_Nest_a_plots[[5]],
          combined_Nest_a_plots[[3]], combined_Nest_a_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1.25,
          ncol = 2, nrow = 3)

invisible(dev.off())
knitr::include_graphics(Nest_a_full_plot)

#Overall model performance 
Baselga_ov <- all_ests$BaselgaMod_Performance
Baselga_ov <- Baselga_ov[,c(1:3, 29:31, 35)] %>% distinct()

Baselga_ov$Estimator <- factor(Baselga_ov$Estimator, levels = baselga_estimators)

Turn_ov <- Baselga_ov[Baselga_ov$Variable == 'Turnover' & (Baselga_ov$Estimator == 'Observed' | Baselga_ov$Type == 'LM'),]
Nest_ov <- Baselga_ov[Baselga_ov$Variable == 'Nestedness' & (Baselga_ov$Estimator == 'Observed' | Baselga_ov$Type == 'LM'),]

Turn_ov_plots <- plotty(Turn_ov,
                        estimators = baselga_estimators,
                        metrics = 'Mean_Accuracy',
                        scenarios = 1:6)

Turn_ov_full_plot <- fs::path('./Figures/',  "Turnover_Model_Overall_Accuracy.png")
agg_png(Turn_ov_full_plot, width = 600, height = 400, units = "mm", res = 500) 

ggarrange(Turn_ov_plots[[1]], Turn_ov_plots[[4]],
          Turn_ov_plots[[2]], Turn_ov_plots[[5]],
          Turn_ov_plots[[3]], Turn_ov_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1,
          ncol = 2, nrow = 3,
          heights = c(1.05, 1, 1))


invisible(dev.off())
knitr::include_graphics(Turn_ov_full_plot)

Nest_ov_plots <- plotty(Nest_ov,
                        estimators = baselga_estimators,
                        metrics = 'Mean_Accuracy',
                        scenarios = 1:6)

Nest_ov_full_plot <- fs::path('./Figures/',  "Nestedness_Model_Overall_Accuracy.png")
agg_png(Nest_ov_full_plot, width = 600, height = 400, units = "mm", res = 500) 

ggarrange(Nest_ov_plots[[1]], Nest_ov_plots[[4]],
          Nest_ov_plots[[2]], Nest_ov_plots[[5]],
          Nest_ov_plots[[3]], Nest_ov_plots[[6]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.5',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 2.0',
                     'Scenario 3: \u03bc = 0.025, \u03c3 = 0.5',
                     'Scenario 4: \u03bc = 0.025, \u03c3 = 2.0',
                     'Scenario 5: \u03bc = 0.1, \u03c3 = 0.5',
                     'Scenario 6: \u03bc = 0.1, \u03c3 = 2.0'), 
          font.label = list(size = 18, face = 'bold.italic'),
          hjust = -1.1, vjust = 1,
          ncol = 2, nrow = 3,
          heights = c(1.05, 1, 1))


invisible(dev.off())
knitr::include_graphics(Nest_ov_full_plot)
