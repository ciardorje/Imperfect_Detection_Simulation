rm(list=ls()); gc()

library(pacman)
p_load(dplyr, plyr)

######Read estimate files and regroup repeat sampling processes######

sim_ests_files <- list.files('./Results', full.names = T)
sim_ests <- vector('list', length(sim_ests_files)/10)
ids <- c()

for(i in 1:length(sim_ests_files)){
 
  load(sim_ests_files[i])
  id <- sim_i_estimates[[1]]$id
  
  if(!id %in% ids){
    
    sim_ests[[id]] <- sim_i_estimates[[1]][1:2]
    sim_ests[[id]]$REsts <- sim_i_estimates[[1]][[4]]
    sim_ests[[id]]$BEsts <- sim_i_estimates[[1]][[5]]
    sim_ests[[id]]$SAREsts <- sim_i_estimates[[1]][[6]]
    sim_ests[[id]]$BModEsts <- sim_i_estimates[[1]][[7]]
    sim_ests[[id]]$BaselgaEsts <- sim_i_estimates[[1]][[8]]
    sim_ests[[id]]$BaselgaModEsts <- sim_i_estimates[[1]][[9]]
    
  } else {
    
    sim_ests[[id]]$REsts <- rbind(sim_ests[[id]]$REsts, sim_i_estimates[[1]][[4]])
    sim_ests[[id]]$BEsts <- rbind(sim_ests[[id]]$BEsts, sim_i_estimates[[1]][[5]])
    sim_ests[[id]]$SAREsts <- rbind(sim_ests[[id]]$SAREsts, sim_i_estimates[[1]][[6]])
    sim_ests[[id]]$BModEsts <- rbind(sim_ests[[id]]$BModEsts, sim_i_estimates[[1]][[7]])
    sim_ests[[id]]$BaselgaEsts <- rbind(sim_ests[[id]]$BaselgaEsts, sim_i_estimates[[1]][[8]])
    sim_ests[[id]]$BaselgaModEsts <- rbind(sim_ests[[id]]$BaselgaModEsts, sim_i_estimates[[1]][[9]])
    
  }
  
  ids <- c(ids, id)
 
}

#-----Evaluate Estimator Performance-----

est_performance <- function(x){
  
  id <- x$id
  scen <- x$Parameters$Detection_Scenario
  params <- x$Parameters
  REsts <- x$REsts
  BEsts <- x$BEsts
  BaselgaEsts <- x$BaselgaEsts
  SAREsts <- x$SAREsts
  BModEsts <- x$BModEsts
  BaselgaModEsts <- x$BaselgaModEsts 
  
  #############Performance in Richness estimates###############
  REsts$Bias <- (REsts$Estimate - REsts$True_R)/REsts$True_R
  REsts$Accuracy <- abs(REsts$Bias)
  REsts$Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                           x = REsts$True_R, 
                                           left = REsts$LCI, 
                                           right = REsts$UCI))
  
  REsts <- REsts %>% 
            group_by(Estimate_ID) %>%
            mutate(Mean_Bias = mean(Bias),
                   Mean_Accuracy = mean(Accuracy),
                   Mean_Coverage = mean(Coverage, na.rm = T),
                   Precision = sd(Estimate)/mean(Estimate),
                   Non_NA_Site = mean(!is.na(SE.SD)),
                   ID = id,
                   Sampling_Events = params$Sampling_Events,
                   bArea_Shape = params$bArea_Shape,
                   mu_ip = params$mu_ip,
                   sd_ip = params$sd_ip,
                   Detection_Scenario = scen)
          
  
  ###############Performance in Beta estimates##################
  BEsts$Bias <- ((BEsts$Estimate - BEsts$True_B) / BEsts$True_B)
  BEsts$Accuracy <- abs(BEsts$Bias)
  BEsts$Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                        x = BEsts$True_B, 
                                        left = BEsts$LCI, 
                                        right = BEsts$UCI))
  
  BEsts <- BEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(Mean_Bias = mean(Bias),
           Mean_Accuracy = mean(Accuracy),
           Mean_Coverage = mean(Coverage, na.rm = T),
           Precision = sd(Estimate)/mean(Estimate),
           Non_NA_Pair = mean(!is.na(SE.SD)),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
  
  ###Baselga estimates###
  BaselgaEsts$Turn_Bias <- ((BaselgaEsts$Turnover - BaselgaEsts$True_Turnover) / BaselgaEsts$True_Turnover)
  BaselgaEsts$Turn_Accuracy <- abs(BaselgaEsts$Turn_Bias)
  BaselgaEsts$Turn_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                      x = BaselgaEsts$True_Turnover, 
                                      left = BaselgaEsts$T_LCI, 
                                      right = BaselgaEsts$T_UCI))
  
  BaselgaEsts$Nest_Bias <- ((BaselgaEsts$Nestedness - BaselgaEsts$True_Nestedness) / BaselgaEsts$True_Nestedness)
  BaselgaEsts$Nest_Accuracy <- abs(BaselgaEsts$Nest_Bias)
  BaselgaEsts$Nest_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                                 x = BaselgaEsts$True_Nestedness, 
                                                 left = BaselgaEsts$N_LCI, 
                                                 right = BaselgaEsts$N_UCI))
  
  BaselgaEsts <- BaselgaEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(Mean_Turn_Bias = mean(Turn_Bias),
           Mean_Turn_Accuracy = mean(Turn_Accuracy),
           Mean_Turn_Coverage = mean(Turn_Coverage, na.rm = T),
           Turn_Precision = sd(Turnover)/mean(Turnover),
           Mean_Nest_Bias = mean(Nest_Bias),
           Mean_Nest_Accuracy = mean(Nest_Accuracy),
           Mean_Nest_Coverage = mean(Nest_Coverage, na.rm = T),
           Nest_Precision = sd(Nestedness)/mean(Nestedness),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
    
  
  #################Performance in SAR Parameters#################
  SAREsts$z_Bias <- ((SAREsts$z - SAREsts$True_z) / SAREsts$True_z)
  SAREsts$z_Accuracy <- abs(SAREsts$z_Bias)
  SAREsts$z_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = SAREsts$True_z, 
                                          left = SAREsts$z_LCI, 
                                          right = SAREsts$z_UCI))
  
  SAREsts$c_Bias <- ((SAREsts$c - SAREsts$True_c) / SAREsts$True_c)
  SAREsts$c_Accuracy <- abs(SAREsts$c_Bias)
  SAREsts$c_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = SAREsts$True_c, 
                                          left = SAREsts$c_LCI, 
                                          right = SAREsts$c_UCI))
  
  SAREsts$Model_Accuracy <- SAREsts$z_Accuracy + SAREsts$c_Accuracy
  
  Incomplete_SARs <- SAREsts[!is.na(SAREsts$Non_NA) & SAREsts$Non_NA < 0.8,] #Remove SARs where richness could only be estimated for less than 80% of patches
  SAREsts <- SAREsts[is.na(SAREsts$Non_NA) | SAREsts$Non_NA >= 0.8,]
  
  SAREsts <- SAREsts %>% 
    group_by(Estimate_ID) %>%
    mutate(z_Mean_Bias = mean(z_Bias),
           z_Mean_Accuracy = mean(z_Accuracy),
           z_Mean_Coverage = mean(z_Coverage, na.rm = T),
           z_Precision = abs(sd(z)/mean(z)),
           c_Mean_Bias = mean(c_Bias),
           c_Mean_Accuracy = mean(c_Accuracy),
           c_Mean_Coverage = mean(c_Coverage, na.rm = T),
           c_Precision = abs(sd(c)/mean(c)),
           Mean_Model_Accuracy = mean(Model_Accuracy),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
  
  SAREsts <- plyr::rbind.fill(SAREsts, Incomplete_SARs)  
  
  SAREsts <- SAREsts %>% 
    group_by(Estimate_ID) %>%
    mutate(Non_NA_Repeats = mean(!is.na(z_Mean_Bias)),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
    
  #################Performance in BetaMod Parameters#################
  
  BModEsts$b_Bias <- ((BModEsts$b - BModEsts$True_b) / BModEsts$True_b)
  BModEsts$b_Accuracy <- abs(BModEsts$b_Bias)
  BModEsts$b_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = BModEsts$True_b, 
                                          left = BModEsts$b_LCI, 
                                          right = BModEsts$b_UCI))
  
  BModEsts$a_Bias <- ((BModEsts$a - BModEsts$True_a) / BModEsts$True_a)
  BModEsts$a_Accuracy <- abs(BModEsts$a_Bias)
  BModEsts$a_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = BModEsts$True_a, 
                                          left = BModEsts$a_LCI, 
                                          right = BModEsts$a_UCI))
  
  BModEsts$Model_Accuracy <- BModEsts$b_Accuracy + BModEsts$a_Accuracy
  
  Incomplete_BMods <- BModEsts[!is.na(BModEsts$Non_NA) & BModEsts$Non_NA < 0.8,]
  BModEsts <- BModEsts[is.na(BModEsts$Non_NA) | BModEsts$Non_NA >= 0.8,]
  
  BModEsts <- BModEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(b_Mean_Bias = mean(b_Bias),
           b_Mean_Accuracy = mean(b_Accuracy),
           b_Mean_Coverage = mean(b_Coverage, na.rm = T),
           b_Precision = abs(sd(b)/mean(b)),
           a_Mean_Bias = mean(a_Bias),
           a_Mean_Accuracy = mean(a_Accuracy),
           a_Mean_Coverage = mean(a_Coverage, na.rm = T),
           a_Precision = abs(sd(a)/mean(a)),
           Mean_Model_Accuracy = mean(Model_Accuracy))
  
  BModEsts <- plyr::rbind.fill(BModEsts, Incomplete_BMods)
    
  BModEsts <- BModEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(Non_NA_Repeats = mean(!is.na(b_Mean_Bias)),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
  
  #Baselga models
  BaselgaModEsts$b_Bias <- ((BaselgaModEsts$b - BaselgaModEsts$True_b) / BaselgaModEsts$True_b)
  BaselgaModEsts$b_Accuracy <- abs(BaselgaModEsts$b_Bias)
  BaselgaModEsts$b_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                                 x = BaselgaModEsts$True_b, 
                                                 left = BaselgaModEsts$b_LCI, 
                                                 right = BaselgaModEsts$b_UCI))
  
  BaselgaModEsts$a_Bias <- ((BaselgaModEsts$a - BaselgaModEsts$True_a) / BaselgaModEsts$True_a)
  BaselgaModEsts$a_Accuracy <- abs(BaselgaModEsts$a_Bias)
  BaselgaModEsts$a_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                                 x = BaselgaModEsts$True_a, 
                                                 left = BaselgaModEsts$a_LCI, 
                                                 right = BaselgaModEsts$a_UCI))
  
  BaselgaModEsts$Model_Accuracy <- BaselgaModEsts$b_Accuracy + BaselgaModEsts$a_Accuracy
  
  BaselgaModEsts <- BaselgaModEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(b_Mean_Bias = mean(b_Bias),
           b_Mean_Accuracy = mean(b_Accuracy),
           b_Mean_Coverage = mean(b_Coverage, na.rm = T),
           b_Precision = abs(sd(b)/mean(b)),
           a_Mean_Bias = mean(a_Bias),
           a_Mean_Accuracy = mean(a_Accuracy),
           a_Mean_Coverage = mean(a_Coverage, na.rm = T),
           a_Precision = abs(sd(a)/mean(a)),
           Mean_Model_Accuracy = mean(Model_Accuracy),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
  
  #################Collate Results######################
  results <- list(Richness_Performance = REsts, 
                  Beta_Performance = BEsts,
                  Baselga_Performance = BaselgaEsts,
                  SAR_Performance = SAREsts,
                  BetaMod_Performance = BModEsts,
                  BaselgaMod_Performance = BaselgaModEsts)
  
  return(results)
  
}


######Apply performance criteria#####
#sim_ests <- sim_ests[lapply(sim_ests, length) > 0]
performance_results <- lapply(sim_ests, function(x) est_performance(x))
all_ests <- do.call(Map, c(f = rbind, performance_results))
save(all_ests, file = 'Estimator_Performance_Results.RData')
