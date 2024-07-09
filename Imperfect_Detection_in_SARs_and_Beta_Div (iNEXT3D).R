rm(list=ls()); gc()

library(pacman)
p_load(dplyr, iNEXT.beta3D, mobsim, nimble, data.table, tidyr, betapart, gridExtra, grid,
       doSNOW, parallel, sf, sp, betafunctions, ggplot2, ggpubr, ragg, ggeffects)

setwd('C:/Users/ciarn/Desktop/PhD/Imperfect Detection Simulation/')

#-----Simulation Parameters-----

n_patches <- 25  #Number of patches in each landscape
max_patch <- 20000  #Maximum patch size in Ha
min_patch <- 25 #Defined smallest patch size to give good range

patch_visits <- c(3, 6, 12)  #Number of visits to each patch

mu_ip <- c(0.005, 0.05, 0.2) #Mean species-individual detection probability
sd_ip <- c(0.25, 1, 3)  #SD of species-individual detection probability
detection_scenarios <- expand.grid(mu_ip, sd_ip)
colnames(detection_scenarios) <- c('Mean_ip', 'SD_ip')
detection_scenarios$Detection_Scenario <- 1:nrow(detection_scenarios) 

#Calculate search area for a 500m transect with 50 search radius
transect_start <- st_as_sf(SpatialPoints(cbind(-56, -9), 
                                         proj4string = CRS('EPSG:4326')))
buffer <- st_buffer(transect_start, 500) #Generate points 500 metres from origin
transect_end <- st_as_sf(SpatialPoints(cbind(buffer[[1]][[1]][[1]][1,1],
                                             buffer[[1]][[1]][[1]][1,2]), 
                                       proj4string = CRS('EPSG:4326')))
transect <- st_cast(st_union(transect_start, transect_end), 'LINESTRING') #Create theoretical transect line

search_area <- st_area(st_buffer(transect, 50))/10000 #Calculate search area in Ha

#N transects sampled in each patch - increases with area
n_transects <- c(rep(1, 5), rep(2, 5), rep(3, 5),
                 rep(4, 5), rep(5, 5))

#Maximum shape parameter for species-level Area response beta distributions
#Higher values yield communities with greater SAR slopes
bArea_Shape <- c(4, 8, 12) 

#Construct all parameter combinations and replicates
params <- expand.grid(patch_visits, bArea_Shape, mu_ip, sd_ip) #Find all parameter combos
params <- params[rep(seq_len(nrow(params)), 35), ] #Replicate param combos 35 times
colnames(params) <- c('Sample_Events', 'bArea_Shape', 'Mean_ip', 'SD_ip')
params <- merge(params, detection_scenarios)

#-----Data Simulation-----
data_sim <- function(patch_visits, bArea_Shape, mu_ip, sd_ip, scen, id){
  
  ####################Simulate patch sizes########################
  #Simulate landscape with many small/medium patches and a few large to represent
  #realistic fragmented habitats. Uses four-parameter beta distribution 
  patches <- ceiling(rBeta.4P(n = (n_patches-2), 
                              l = min_patch, u = max_patch,  
                              alpha = 1, beta = 4))
  patches <- c(min_patch, patches, max_patch) #defined min (10Ha) and max (20000Ha) patches
  patches <- patches[order(patches)]
  #hist(patches)
  
  
  ##################Simulate landscape metacommunity#################
  
  ###Simulate landscape-level species abundance distribution (log-normal)###
  density <- 500 #Number of organisms / Ha
  nSpp <- 200 #No. of simulated spp
  indis <- sum(patches) * density  #Total pop size = total patch area * density
  sad <- sim_sad(s_pool = nSpp, n_sim = indis/density, sad_type = 'lnorm', 
                 sad_coef = list('meanlog' = 650, 'sdlog' = 3), fix_s_sim = T) * density
  #Generating SAD for capacity/density then * by density ensures a viable minimum
  #landscape population for each species (i.e., 500)
  #sad
  
  ###Assign individuals to patches###
  pop_sizes <- (patches * density) #Determine pop size of each patch according to area
  patch_pops <- data.frame(Patch = 1:n_patches, 
                           Area = patches, 
                           Population = pop_sizes)
  
  #Convert SAD to species x count matrix (landscape-level)
  spp_pot <- data.frame(Species = names(sad),
                        Count = as.numeric(sad)) 
  spp_pot <- spp_pot[order(spp_pot$Count, decreasing = T),]
  spp_pot$Species <- 1:nSpp
  #hist(spp_pot$Count)
  
  #Assign species with weighted probabilities according to generated area responses
  #Generate positive and negative responses with four-parameter beta distribution
  spp_pot$alpha <- seq(0.1, bArea_Shape, length.out = nSpp) #Species level shape parameters, 
  #increasing chance of severe area response with increasing rarity
  #ensures realistic, viable patch-level pops of rarest species where they occur
  spp_pot$bArea <- round(mapply(FUN = rBeta.4P,                              
                                n = 1,
                                u = 5, 
                                l = 0, 
                                alpha = spp_pot$alpha,
                                beta = 5), 3) 
  
  #Select 1/8th of species to have negative area responses
  nNegSpp <- round(nSpp/8)
  negSpp <- sample(1:nSpp, nNegSpp) #Randomly select 1/6th of species
  spp_pot$bArea[negSpp] <- -spp_pot$bArea[negSpp] #Assign generated negative response
  
  #Generate probabilities of assignment to each patch for spp with defined area response
  assignment_probs <- matrix(NA, nrow = nSpp, ncol = n_patches)
  
  for(i in 1:nSpp){
    for(x in 1:n_patches){
      #Probability of assignment to patch = Exponent of area response * log(Patch Area)
      assignment_probs[i,x] <- exp(spp_pot$bArea[i] * log(patches[x]))
    }
    #Scale to sum to 1
    assignment_probs[i,] <- (assignment_probs[i,] / sum(assignment_probs[i,])) 
  }
  
  #Create Site x Species Matrix
  true_abun <- matrix(0, nrow = n_patches, ncol = nSpp)
  
  locs <- 1:n_patches #List of fragment IDs (used to remove fragments from pool once capacity reached)
  done <- FALSE
  
  while(sum(spp_pot$Count) > 0){ #While individuals in the species pool remain to be assigned...
    for(sp in 1:nSpp){ #Repeatedly assign one individual from each species until all done
      if(spp_pot$Count[sp] > 0){ #If species individuals remain to be assigned proceed
        vacancies <- 0 #Set vacancies to 0 to enter next loop
        while(vacancies == 0){ #If there are 0 vacancies in patch selected in next section, retry assignment
          patch <- sample(locs, 1, prob = assignment_probs[sp, locs]) #Randomly select patch according to generated assignment probabilities
          vacancies <- patch_pops$Population[patch] #Find no. of individual vacancies yet to be filled in the patch
          if(vacancies > 0){ #If vacancies available in patch...
            true_abun[patch, sp] <- (true_abun[patch, sp] + 100) #Assign an individual to the patch
            spp_pot$Count[sp] <- (spp_pot$Count[sp] - 100) #Decrease species count by 1
            patch_pops$Population[patch] <- (patch_pops$Population[patch] - 100) #Decrease remaining vacancies in patch by 1
          } else { #If patch is already full...
            locs <- locs[!locs %in% patch] #Remove patch from selection process 
          }
          if(length(locs) == 1){ #If only one patch has vacancies...
            for(i in 1:nrow(spp_pot)){ #Assign all remaining individuals to that patch
              true_abun[locs, i] <- (true_abun[locs, i] + spp_pot$Count[i])
            }
            spp_pot$Count <- 0 #Set population to 0 (all individuals have been assigned)
            done <- TRUE #Convert 'done' to true to end loop
            break
          }
        }
      }
    }
    if(sum(spp_pot$Count) %% 10000 == 0){cat('\n', sum(spp_pot$Count))}
    if(done){break}
  }
  
  #Determine true richness at each site
  true_pa <- ifelse(true_abun > 0, 1, 0) #Presence absence matrix
  True_R <- rowSums(true_pa, na.rm = T)  #Summate spp richness
  patchR <- data.frame(Patch = 1:n_patches, Area = patches, True_R = True_R)
  # plot(patchR$Area, patchR$True_R, ylim = c(0, nSpp))
  # summary(lm(log(True_R) ~ log(Area), data = patchR))
  # 
  #Assess whether species density increases/decreases linearly or curvilinearly with area
  # density <- true_abun
  # for(i in 1:n_patches){
  #   density[i,] <- density[i,]/patches[i]
  # }
  # 
  # for(i in 1:nSpp){
  #   par(mfrow = c(1,2))
  #   plot(patches, true_abun[,i])
  #   plot(patches, density[,i])
  #  }

  
  #####################Simulate sampling#######################
  
  #Incomplete spatial sampling - isolate indis 'present' within sample area
  patch_proportion <- as.numeric(((search_area/patches) * n_transects))
  
  #Species-specific, individual-level detection probabilities
  ip <- plogis(rnorm(nSpp, qlogis(mu_ip), sd_ip))  #(mu & sd on logit scale; back to prob)
  
  observed_data <- vector('list', 10)
  
  #Simulate 10 repeat sampling processes
  for(i in 1:10){
    
    set.seed(i)
    
    #Number of individuals within sampled area
    sampleable <- matrix(NA, nrow = n_patches, ncol = nSpp)
    for(j in 1:n_patches){
      for(x in 1:nSpp){
        #sampleable[j,x] <- rbinom(1, true_abun[j,x], patch_proportion[j]) 
        sampleable[j,x] <- round(true_abun[j,x] * patch_proportion[j])
      }
    }
    
    Sampleable <- rowSums(sampleable > 0)
    #plot(patchR$Area, Sampleable,  ylim = c(0, nSpp))
    #lm(log(Sampleable) ~ log(patchR$Area))
    
    #Sampling visits
    Y <- array(NA, dim = c(n_patches, nSpp, patch_visits))
    for(j in 1:n_patches){
      for(x in 1:nSpp){
        Y[j,x, ] <- rbinom(n = patch_visits, 
                           size = sampleable[j, x], 
                           prob = ip[x])
      }
    }
    
    #Remove unobserved species
    sp_obs <- data.frame(Species = 1:nSpp, Observations = apply(Y, 2, sum))
    observed <- sp_obs$Species[sp_obs$Observations > 0]
    Y <- Y[, observed,]
    
    #Summarise Observed data
    Ypooled <- apply(Y, c(1,2), sum) #Sum species x site abundance across all sample events
    Yinc <- ifelse(Y > 0, 1, 0) #Convert abundance samples to incidence data
    Yinc <- apply(Yinc, c(1,2), sum)
    
    Z <- ifelse(Yinc > 0, 1, 0) #Observed presence absence across all samples
    Observed_R <- rowSums(Z > 0) #Calculate observed patch level richness
    ObsGamma <- length(observed)
    #plot(patchR$Area, Observed_R,  ylim = c(0, nSpp))
    #lm(log(Observed_R) ~ log(patchR$Area))
    
    observed_data[[i]] <- list(Y = Y, Ypooled = Ypooled, Yinc = Yinc, Z = Z, 
                               Observed_R = Observed_R, Observed_Gamma = ObsGamma, Repeat = i)
    
  }
  
  
  
  #########################Save results#############################
  parameters <- list(Sampling_Events = patch_visits, bArea_Shape = bArea_Shape, 
                     mu_ip = mu_ip, sd_ip = sd_ip, Detection_Scenario = scen,
                     patch_prop = patch_proportion)
  
  true_data <- list(Patches = patchR, True_Abundance = true_abun, True_PA = true_pa,
                    Spp_ip = ip, Gamma = nSpp)
  
  data <- list(id = id, Parameters = parameters, 
               True_Data = true_data, Observed_Data = observed_data)
  return(data)
  
}

######Simulate 30 datasets for each parameter combination#####
cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)
ntasks <- nrow(params)
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(x) setTxtProgressBar(pb, x)
opts <- list(progress = progress)

landscape_simulations <- foreach(i = 1:nrow(params), 
                                 .packages = c('mobsim', 
                                               'betafunctions'),
                                 .options.snow = opts,
                                 .inorder = F
  ) %dopar% {
  
  set.seed(i)
  data_sim(patch_visits = params$Sample_Events[i],
           bArea_Shape = params$bArea_Shape[i],
           mu_ip = params$Mean_ip[i], 
           sd_ip = params$SD_ip[i],
           id = i,
           scen = params$Detection_Scenario[i])
  }

stopCluster(cl)

save(landscape_simulations, file = 'Landscape_Community_Sims.RData')
load('Landscape_Community_Sims.RData')

######Parse simulation repeats into individual sampling processes for running on cluster#####
landscape_sims_parsed <- vector('list', (nrow(params)*10))
x <- 1

for(i in 1:length(landscape_simulations)){
  for(n in 1:10){
    
    landscape_sims_parsed[[x]] <- list(id = landscape_simulations[[i]][[1]],
                                       Parameters = landscape_simulations[[i]][[2]],
                                       True_Data = landscape_simulations[[i]][[3]],
                                       Observed_Data = list(
                                       landscape_simulations[[i]][[4]][[n]]))
    x <- x+1
  }
}

save(landscape_sims_parsed, file = 'Landscape_Sims_Parsed.RData')

#-----Apply Richness and Beta Estimators-----

load('Landscape_Sims_Parsed.RData')

#####Estimation Functions#####

#Formula for sorensen index
sorensen <- function(x, y){ ((2 * sum(x * y)) / (sum(x) + sum(y))) }


#Hierarchical Multi-Species Occupancy Model structure
HMSOM <- nimbleCode({
  
  #Community-Level Hyperparameters#
  omega ~ dbeta(0.001, 1) #Community inclusion parameter - approximation of Link's scale prior
  
  psi.mean ~ dbeta(1, 1) #Approximation of Uniform distribution
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  mu.Area ~ dnorm(0,0.1) 
  sd.Area ~ dunif(0, 5)
  tau.Area <- 1/sd.Area^2
  
  p.mean ~ dbeta(1, 1) #Approximation of Uniform distribution
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  rho ~ dunif(-1, 1) 
  tau.eta <- tau.lp/(1 - rho^2)
  
  for(k in 1:M){
    
    #Species-Level Priors#
    w[k] ~ dbern(omega) 
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    bArea[k] ~ dnorm(mu.Area, tau.Area)
    mu.eta[k] <- mu.lp + rho * sd.lp/sd.lpsi * (lpsi[k] - mu.lpsi)
    lp[k] ~ dnorm(mu.eta[k], tau.eta)
    logit(p[k]) <- lp[k]
    
    for(i in 1:nPatches){
      
      #Likelihood#
      logit(psi[i, k]) <- lpsi[k] + bArea[k] * Area[i] #Non-conditional occurrence probability 
      z[i, k] ~ dbern(psi[i, k] * w[k] ) #Binary occurrence indicator
      y[i, k] ~ dbin(p[k] * z[i, k], nSamples)
      
    }
  }
  
  #Derived Richness Estimates
  for(i in 1:nPatches){
    PatchR[i] <- sum(z[i,1:M])
  } 
  
})


####Function to derive richness and beta estimates from simulated data####
estimate_diversity <- function(x){
  
  ###Set-up###
  id <- x[[1]]
  inputs <- x[[2]]
  true_data <- x[[3]]
  obs_data <- x[[4]]
  n_patches <- nrow(true_data$Patches)
  
  results <- vector('list', length(obs_data))
  
  for(i in 1:length(obs_data)){
    
    data <- obs_data[[i]]
    
    #################### Richness #########################
    
    ###Chao###
    pooled_abun <- data$Ypooled #All sampling events pooled
    chao_results <- data.frame(matrix(ncol = 6, nrow = n_patches))
    colnames(chao_results) <- c('Estimate', 'SE.SD', 'LCI', 'UCI', 'Patch', 'Estimator')
    
    for(j in 1:n_patches){
      
      print(j)
      
      #Asymptotic chao richness estimate
      #Doesn't fail with no singletons or doubletons, just returns observed R w/ 0 uncertainty
      result <- AO3D(pooled_abun[j,], q = 0, datatype = 'abundance')[1,3:6] 
      
      result$Patch <- j
      result$Estimator <- 'Chao'
      chao_results[j,] <- result[1,]
      
    }
    
    
    ###HMSOM###
    #Data preparation
    nSamples <- inputs$Sampling_Events
    Area <- true_data$Patches$Area
    
    nSpObs <- data$Observed_Gamma  
    AugSp <- 2 * nSpObs #Augment w/ 2 * no. observed sps
    M <- nSpObs + AugSp
    w <- c(rep(1, nSpObs), rep(NA, AugSp))
    
    y <- data$Yinc 
    yAug <- cbind(y, matrix(0, nrow = n_patches, ncol = AugSp)) #Augment detection array
    z <- ifelse(yAug > 0, 1, NA)
    
    HMSOM_data <- list(y = yAug, z = z, M = M, w = w,
                       nPatches = n_patches, nSamples = nSamples, 
                       Area = as.numeric(scale(Area)))
    
    #Model Fitting
    model <- nimbleModel(HMSOM, constants = HMSOM_data)
    MCMCconf <- configureMCMC(model,
                              monitors = c('PatchR', 'z'))
    MCMC <- buildMCMC(MCMCconf)
    compModel <- compileNimble(model)
    compMCMC <- compileNimble(MCMC, project = compModel)
    HMSOM_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                         thin = 20, nchains = 4, samplesAsCodaMCMC = T)
    
    ##Extract Posterior richness estimates##
    HMSOM_df <- as.data.frame(as.matrix(HMSOM_out))
    HMSOM_R <- HMSOM_df[,grepl('PatchR', names(HMSOM_df))]
    
    #Richness in sampled areas
    HMSOM_R_Mn <- apply(HMSOM_R, 2, mean)
    HMSOM_R_sd <- apply(HMSOM_R, 2, sd)
    HMSOM_R_lci <- apply(HMSOM_R, 2, quantile, probs = 0.025)
    HMSOM_R_uci <- apply(HMSOM_R, 2, quantile, probs = 0.975)
    
    #Collate HMSOM estimates
    HMSOM_R_results <- as.data.frame(cbind(HMSOM_R_Mn, HMSOM_R_sd, HMSOM_R_lci, HMSOM_R_uci))
    HMSOM_R_results$Patch <- 1:n_patches
    HMSOM_R_results$Estimator <- c(rep('HMSOM', n_patches))
    HMSOM_R_results$Samples <- NA
    colnames(HMSOM_R_results)[1:4] <- c('Estimate', 'SE.SD', 'LCI', 'UCI')
    
    ###Observed###
    observed_results <- data.frame(Estimate = rowSums(data$Z > 0),
                                   SE.SD = NA,
                                   LCI = NA,
                                   UCI = NA,
                                   Patch = 1:n_patches, 
                                   Estimator = 'Observed')
    
    #Combine Results
    patchR <- rbind(chao_results, HMSOM_R_results, observed_results)
    patchR <- merge(true_data$Patches, patchR)
    patchR$Estimate_ID <- 1:nrow(patchR)
    patchR$Repeat <- data$Repeat
    
    
    ##################### Beta Diversity ##########################
    
    #Create Pairwise Beta matrices
    
    pairs <- n_patches * (n_patches - 1)/2 #n pairs
    P1 <- P2 <- c()
    
    for(p1 in 1:(n_patches - 1)){ #Create all pairs
      for(p2 in (p1+1):n_patches){
        P1 <- c(P1, p1)
        P2 <- c(P2, p2)
      }
    }
    
    Observed_B <-   
      HMSOM_B_Ests <- 
      Chao_B_Ests <- data.frame(Site_1 = P1,
                                Site_2 = P2,
                                Estimate = NA,
                                SE.SD = NA,
                                LCI = NA,
                                UCI = NA)
    
    #True Beta
    true_comm <- true_data$True_PA
    True_B <- data.frame(Site_1 = P1,
                         Site_2 = P2,
                         True_B = NA)
    
    #Extract HMSOM z matrix
    HMSOM_z <- HMSOM_df[,grepl('z', names(HMSOM_df))]
    
    indices <- vapply(strsplit(colnames(HMSOM_z), "[", fixed = T), 
                      `[`, 2, FUN.VALUE = character(1))
    indices <- vapply(strsplit(indices, "]", fixed = T), 
                      `[`, 1, FUN.VALUE = character(1))
    indices <- strsplit(indices, ",")
    ri <- as.numeric(sapply(indices, '[[', 1))
    ci <- as.numeric(sapply(indices, '[[', 2))
    z_mat <- array(NA, dim = c(max(ri), max(ci), nrow(HMSOM_z)))
    
    for(ind in 1:length(indices)){ 
      z_mat[ri[ind],ci[ind],] <- HMSOM_z[,ind]
    }
    
    #Estimate beta
    for(k in 1:pairs){
      
      #Observed Beta
      site_pair <- as.data.frame(t(pooled_abun[c(P1[k],P2[k]),]))
      Observed_B[k,'Estimate'] <- sorensen(site_pair[,1] > 0, site_pair[,2] > 0)
      
      #True beta
      True_B[k,'True_B'] <- sorensen(true_comm[P1[k],], true_comm[P2[k],])
      
      #Chao Asymptotic Beta (Sorensen Dissim)
      #doesn't fail with no shared sp, just returns 1 w/ 0 uncertainty
      #Estimate asymptotic sorensen dissim
      chao_dissim <- iNEXTbeta3D(site_pair, q = 0, datatype = 'abundance', level = 1)[[1]][['1-C']][,c(5,7:9)]
      #Convert to similarity
      Chao_B_Ests$Estimate[k] <- 1 - chao_dissim[,1]
      Chao_B_Ests$SE.SD[k] <- chao_dissim[,2]
      Chao_B_Ests$LCI[k] <- 1 - chao_dissim[,4]
      Chao_B_Ests$UCI[k] <- 1 - chao_dissim[,3]
      
      #HMSOM Beta 
      sor_ests <- vector('numeric', dim(z_mat)[3])
      
      for(iter in 1:dim(z_mat)[3]){ 
        sor_ests[iter] <- sorensen(z_mat[P1[k],,iter], z_mat[P2[k],,iter]) 
      }
      
      #HMSOM sorensen estimate summary stats
      HMSOM_B_Ests$Estimate[k] <- mean(sor_ests)
      HMSOM_B_Ests$SE.SD[k] <- sd(sor_ests)
      HMSOM_B_Ests$LCI[k] <- quantile(sor_ests, 0.025)
      HMSOM_B_Ests$UCI[k] <- quantile(sor_ests, 0.975)
      
      cat(k, '\n')
    }
    
    HMSOM_B_Ests$Estimator <- 'HMSOM'
    Observed_B$Estimator <- 'Observed'
    Chao_B_Ests$Estimator <- 'Chao'
    
    patchB <- rbind(Observed_B, HMSOM_B_Ests, NPE_B_Ests_Single)
    patchB <- merge(patchB, True_B)
    patchB$Estimate_ID <- 1:nrow(patchB)
    patchB$Repeat <- data$Repeat
    
    
    ######################## SARs ############################
    
    #True SAR
    True_SAR_Data <- unique(patchR[, c('Area', 'True_R')])
    True_SAR <- lm(log(True_R + 1) ~ log(Area), data = True_SAR_Data)
    True_SAR_Summary <- c(True_SAR$coefficients[1], 
                          True_SAR$coefficients[2])
    
    #Observed SAR
    Observed_SAR_Data <- patchR[patchR$Estimator == 'Observed', c('Area', 'Estimate')]
    Observed_SAR <- lm(log(Estimate + 1) ~ log(Area), data = Observed_SAR_Data)
    Observed_SAR_Summary <- c('Observed', 'LM', Observed_SAR$coefficients[1], 
                              confint(Observed_SAR, '(Intercept)', level = 0.95)[1],
                              confint(Observed_SAR, '(Intercept)', level = 0.95)[2],
                              summary(Observed_SAR)[['coefficients']][1,2],
                              Observed_SAR$coefficients[2],
                              confint(Observed_SAR, 'log(Area)', level = 0.95)[1],
                              confint(Observed_SAR, 'log(Area)', level = 0.95)[2],
                              summary(Observed_SAR)[['coefficients']][2,2])
    
    #Estimated SARs
    
    HMSOM_SAR_data <- patchR[patchR$Estimator == 'HMSOM', c('Area', 'Estimate')]
    Chao_SAR_Data <- patchR[patchR$Estimator == 'Chao', c('Area', 'Estimate')]
    
    SARs <- list(
      
      HMSOM = lm(log(Estimate + 1) ~ log(Area), data = HMSOM_SAR_data),
      Chao = lm(log(Estimate + 1) ~ log(Area), data = Chao_SAR_Data),
      
    )
    
    #Collate results
    SAR_results <- as.data.frame(matrix(nrow = 2, ncol = 9))
    colnames(SAR_results) <- c('Estimator', 'c',  'c_LCI', 'c_UCI', 'c_SE',
                               'z', 'z_LCI', 'z_UCI', 'z_SE')
    
    for(n in 1:length(SARs)){
      
      SAR_results$Estimator[n] <- names(SARs)[n]
      SAR_results$c[n] <- as.numeric(coef(SARs[[n]])[1])
      SAR_results$z[n] <- as.numeric(coef(SARs[[n]])[2])
      SAR_results$c_LCI[n] <- confint(SARs[[n]], '(Intercept)', 0.95)[1]
      SAR_results$c_UCI[n] <- confint(SARs[[n]], '(Intercept)', 0.95)[2]
      SAR_results$c_SE[n] <- summary(SARs[[n]])[['coefficients']][1,2]
      SAR_results$z_LCI[n] <- confint(SARs[[n]], 'log(Area)', 0.95)[1]
      SAR_results$z_UCI[n] <- confint(SARs[[n]], 'log(Area)', 0.95)[2]
      SAR_results$z_SE[n] <- summary(SARs[[n]])[['coefficients']][2,2]
      
    }
    
    SAR_results <- rbind(Observed_SAR_Summary, SAR_results)
    SAR_results$True_z <-  True_SAR_Summary[2]
    SAR_results$True_c <- True_SAR_Summary[1]
    SAR_results$Estimate_ID <- 1:nrow(SAR_results)
    SAR_results$Repeat <- data$Repeat
    SAR_results[,3:13] <- apply(SAR_results[,3:14], 2, as.numeric)
    
    
    ################ Beta Diversity Models ###################
    
    #Create environmental distance matrix
    dist_mat <- matrix(NA, ncol = n_patches, nrow = n_patches)
    
    for(n in 1:n_patches){
      for(o in 1:n_patches){
        if(n != o){
          
          dist_mat[n,o] <- abs(true_data$Patches$Area[n] - true_data$Patches$Area[o])
          
        } 
      }
    }
    
    dist_mat <- log(dist_mat + 1)
    dist_mat <- reshape2::melt(dist_mat, varnames = c('Site_2', 'Site_1'))
    dist_mat <- dist_mat[dist_mat$Site_1 != dist_mat$Site_2,]
    dist_mat <- dist_mat[!duplicated(t(apply(dist_mat[c("Site_2", "Site_1")], 1, sort))), ]
    colnames(dist_mat)[3] <- 'Area_Diff'
    
    
    #Merge with beta estimates
    BModDat <- merge(patchB, dist_mat)
    BModDat$Estimate[BModDat$Estimate == 1 | BModDat$Estimate > 1] <- 1-1e-3
    BModDat$Estimate[BModDat$Estimate == 0 | BModDat$Estimate < 0] <- 0+1e-3
    BModDat$True_B[BModDat$True_B == 1] <- 1-1e-3
    BModDat$True_B[BModDat$True_B == 0] <- 0+1e-3
    
    #True Beta model
    Btrue <- unique(BModDat[,c('Area_Diff', 'True_B')])
    Btrue_mod <- glm(True_B ~ Area_Diff, family = binomial, data = Btrue)
    Btrue_mod_summ <- c(Btrue_mod$coefficients[1], 
                        Btrue_mod$coefficients[2])
    
    #Observed Beta model
    Bobs <- BModDat[BModDat$Estimator == 'Observed', c('Area_Diff', 'Estimate')]
    Bobs_mod <- glm(Estimate ~ Area_Diff, family = binomial, data = Bobs)
    Bobs_mod_summ <- c('Observed', 'LM', Bobs_mod$coefficients[1], 
                       confint(Bobs_mod, '(Intercept)', level = 0.95)[1],
                       confint(Bobs_mod, '(Intercept)', level = 0.95)[2],
                       summary(Bobs_mod)[['coefficients']][1,2],
                       Bobs_mod$coefficients[2],
                       confint(Bobs_mod, 'Area_Diff', level = 0.95)[1],
                       confint(Bobs_mod, 'Area_Diff', level = 0.95)[2],
                       summary(Bobs_mod)[['coefficients']][2,2])
    
    #Estimator Beta models
    Bchao_dat <- BModDat[BModDat$Estimator == 'Chao', c('Area_Diff', 'Estimate')]
    Bhmsom_dat <- BModDat[BModDat$Estimator == 'HMSOM', c('Area_Diff', 'Estimate')]
    
    Bmods <- list(
      
      Chao = glm(Estimate ~ Area_Diff, family = binomial, data = Bchao_dat),
      HMSOM = glm(Estimate ~ Area_Diff, family = binomial, data = Bhmsom_dat)
      
    )
    
    #Collate results
    BMod_Results <- as.data.frame(matrix(nrow = 3, ncol = 9))
    colnames(BMod_Results) <- c('Estimator', 'a',  'a_LCI', 'a_UCI', 'a_SE',
                                'b', 'b_LCI', 'b_UCI', 'b_SE')
    
    for(n in 1:length(Bmods)){
      
      BMod_Results$Estimator[n] <- sub('x', '', names(Bmods)[n])
      BMod_Results$a[n] <- as.numeric(coef(Bmods[[n]])[1])
      BMod_Results$b[n] <- as.numeric(coef(Bmods[[n]])[2])
      BMod_Results$a_LCI[n] <- confint(Bmods[[n]], '(Intercept)', 0.95)[1]
      BMod_Results$a_UCI[n] <- confint(Bmods[[n]], '(Intercept)', 0.95)[2]
      BMod_Results$a_SE[n] <- summary(Bmods[[n]])[['coefficients']][1,2]
      BMod_Results$b_LCI[n] <- confint(Bmods[[n]], 'Area_Diff', 0.95)[1]
      BMod_Results$b_UCI[n] <- confint(Bmods[[n]], 'Area_Diff', 0.95)[2]
      BMod_Results$b_SE[n] <- summary(Bmods[[n]])[['coefficients']][2,2]
      
    }
    
    BMod_Results <- rbind(BMod_Results, Bobs_mod_summ)
    BMod_Results$True_a <- Btrue_mod_summ[1]
    BMod_Results$True_b <- Btrue_mod_summ[2]
    BMod_Results$Estimate_ID <- 1:nrow(BMod_Results)
    BMod_Results$Repeat <- data$Repeat
    BMod_Results[,3:13] <- apply(BMod_Results[,3:13], 2, as.numeric)
    
    #####Collate all results and save#####
    results[[i]] <- list(id = id, Parameters = inputs, Repeat = data$Repeat, 
                         Richness_Estimates = patchR, Beta_Estimates = patchB,
                         SAR_Estimates = SAR_results, BetaMod_Estimates = BMod_Results)
    
    
  }
  
  return(results)
  
}


#####Apply estimators to simulated data#####

#This section was not run and a similar process was instead run on cluster
#Below code is just to give an idea of the output structure

#Parallel processing
cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)
ntasks <- length(landscape_sims_parsed)
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(x) setTxtProgressBar(pb, x)
opts <- list(progress = progress)

sim_ests <- foreach(i = 1:length(landscape_sims_parsed),
                     .options.snow = opts,
                     .inorder = F,
                     .packages = c('nimble', 'iNEXT.3D', 'iNEXT.beta3D', 'reshape2')
  ) %dopar% {
  
  out <- estimate_diversity(landscape_sims_parsed[i])
  return(out)
  
}

stopCluster(cl)

for(i in 1:length(sim_ests)){
  save(sim_ests[[i]], file = paste0('./Results/Comm_Ests_', i, '.RData'))
}

gc()


######Read estimate files and regroup repeat sampling processes######

sim_ests_files <- list.files('./Results', full.names = T)
sim_ests <- vector('list', ceiling(length(sim_ests_files)/10))
ids <- c()
n <- 1

for(i in 1:length(sim_ests_files)){
 
  load(sim_ests_files[i])
  id <- as.character(out$id)
  
  if(!id %in% ids){
    
    names(sim_ests)[n] <- id
    n <- n + 1
    ids <- c(ids, id)
    
    sim_ests[[id]] <- out[1:2]
    sim_ests[[id]]$REsts <- out[[4]]
    sim_ests[[id]]$BEsts <- out[[5]]
    sim_ests[[id]]$SAREsts <- out[[6]]
    sim_ests[[id]]$BModEsts <- out[[7]]
    
  } else {
    
    sim_ests[[id]]$REsts <- rbind(sim_ests[[id]]$REsts, out[[4]])
    sim_ests[[id]]$BEsts <- rbind(sim_ests[[id]]$BEsts, out[[5]])
    sim_ests[[id]]$SAREsts <- rbind(sim_ests[[id]]$SAREsts, out[[6]])
    sim_ests[[id]]$BModEsts <- rbind(sim_ests[[id]]$BModEsts, out[[7]])
    
  }
  
  print(i)
 
}

save(sim_ests, file = 'Regrouped_Estimates.RData')

#-----Evaluate Estimator Performance-----

load('Regrouped_Estimates.RData')

est_performance <- function(x){
  
  id <- x$id
  scen <- x$Parameters$Detection_Scenario
  params <- x$Parameters
  REsts <- x$REsts
  BEsts <- x$BEsts
  SAREsts <- x$SAREsts
  BModEsts <- x$BModEsts
  
  #############Performance in Richness estimates###############
  REsts$Bias <- (REsts$Estimate - REsts$True_R)/REsts$True_R
  REsts$Accuracy <- abs(REsts$Bias)
  REsts$Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                           x = REsts$True_R, 
                                           left = REsts$LCI, 
                                           right = REsts$UCI))
  REsts$Width95 <- (REsts$UCI - REsts$LCI) / REsts$Estimate
  
  REsts <- REsts %>% 
            group_by(Estimate_ID) %>%
            mutate(Mean_Bias = mean(Bias),
                   Mean_Accuracy = mean(Accuracy) * -1,
                   Mean_Coverage = mean(Coverage, na.rm = T),
                   Precision = (sd(Estimate)/mean(Estimate)) * -1,
                   Mean_Width_95 = mean(Width95, na.rm = T),
                   ID = id,
                   Sampling_Events = params$Sampling_Events,
                   bArea_Shape = params$bArea_Shape,
                   mu_ip = params$mu_ip,
                   sd_ip = params$sd_ip,
                   Detection_Scenario = scen)
          
  
  ###############Performance in Beta estimates##################
  BEsts$Bias <- BEsts$Estimate - BEsts$True_B
  BEsts$Accuracy <- abs(BEsts$Bias)
  BEsts$Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                        x = BEsts$True_B, 
                                        left = BEsts$LCI, 
                                        right = BEsts$UCI))
  BEsts$Width95 <- BEsts$UCI - BEsts$LCI
  
  BEsts <- BEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(Mean_Bias = mean(Bias),
           Mean_Accuracy = mean(Accuracy) * -1,
           Mean_Coverage = mean(Coverage, na.rm = T),
           Precision = (sd(Estimate)/mean(Estimate)) * -1,
           Mean_Width95 = mean(Width95, na.rm = T),
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
    
  
  #################Performance in SAR Parameters#################
  SAREsts$z_Bias <- SAREsts$z - SAREsts$True_z
  SAREsts$z_Accuracy <- abs(SAREsts$z_Bias)
  SAREsts$z_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = SAREsts$True_z, 
                                          left = SAREsts$z_LCI, 
                                          right = SAREsts$z_UCI))
  SAREsts$z_Width95 <- SAREsts$z_UCI - SAREsts$z_LCI
  
  SAREsts$c_Bias <- SAREsts$c - SAREsts$True_c
  SAREsts$c_Accuracy <- abs(SAREsts$c_Bias)
  SAREsts$c_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = SAREsts$True_c, 
                                          left = SAREsts$c_LCI, 
                                          right = SAREsts$c_UCI))
  SAREsts$c_Width95 <- SAREsts$c_UCI - SAREsts$c_LCI
  
  SAREsts$Model_Accuracy <- SAREsts$z_Accuracy + SAREsts$c_Accuracy
  
  SAREsts <- SAREsts %>% 
    group_by(Estimate_ID) %>%
    mutate(z_Mean_Bias = mean(z_Bias),
           z_Mean_Accuracy = mean(z_Accuracy) * -1,
           z_Mean_Coverage = mean(z_Coverage, na.rm = T),
           z_Precision = (abs(sd(z)/mean(z))) * -1,
           z_Mean_Width95 = mean(z_Width95, na.rm = T),
           c_Mean_Bias = mean(c_Bias),
           c_Mean_Accuracy = mean(c_Accuracy) * -1,
           c_Mean_Coverage = mean(c_Coverage, na.rm = T),
           c_Precision = (abs(sd(c)/mean(c))) * -1,
           c_Mean_Width95 = mean(c_Width95, na.rm = T),
           Mean_Model_Accuracy = mean(Model_Accuracy) * -1,
           ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
  
  SAREsts <- SAREsts %>% 
    group_by(Estimate_ID) %>%
    mutate(ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
    
  #################Performance in BetaMod Parameters#################
  
  BModEsts$b_Bias <- BModEsts$b - BModEsts$True_b
  BModEsts$b_Accuracy <- abs(BModEsts$b_Bias)
  BModEsts$b_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = BModEsts$True_b, 
                                          left = BModEsts$b_LCI, 
                                          right = BModEsts$b_UCI))
  BModEsts$b_Width95 <- BModEsts$b_UCI - BModEsts$b_LCI
  
  BModEsts$a_Bias <- BModEsts$a - BModEsts$True_a
  BModEsts$a_Accuracy <- abs(BModEsts$a_Bias)
  BModEsts$a_Coverage <- as.integer(mapply(FUN = dplyr::between, 
                                          x = BModEsts$True_a, 
                                          left = BModEsts$a_LCI, 
                                          right = BModEsts$a_UCI))
  BModEsts$a_Width95 <- BModEsts$a_UCI - BModEsts$a_LCI
  
  BModEsts$Model_Accuracy <- BModEsts$b_Accuracy + BModEsts$a_Accuracy
  
  BModEsts <- BModEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(b_Mean_Bias = mean(b_Bias),
           b_Mean_Accuracy = mean(b_Accuracy) * -1,
           b_Mean_Coverage = mean(b_Coverage, na.rm = T),
           b_Precision = (abs(sd(b)/mean(b))) * -1,
           b_Mean_Width95 = mean(b_Width95, na.rm = T),
           a_Mean_Bias = mean(a_Bias),
           a_Mean_Accuracy = mean(a_Accuracy) * -1,
           a_Mean_Coverage = mean(a_Coverage, na.rm = T),
           a_Precision = (abs(sd(a)/mean(a))) * -1,
           a_Mean_Width95 = mean(a_Width95, na.rm = T),
           Mean_Model_Accuracy = mean(Model_Accuracy) * -1)
    
  BModEsts <- BModEsts %>% 
    group_by(Estimate_ID) %>%
    mutate(ID = id,
           Sampling_Events = params$Sampling_Events,
           bArea_Shape = params$bArea_Shape,
           mu_ip = params$mu_ip,
           sd_ip = params$sd_ip,
           Detection_Scenario = scen)
  
  
  #################Collate Results######################
  results <- list(Richness_Performance = REsts, 
                  Beta_Performance = BEsts,
                  SAR_Performance = SAREsts,
                  BetaMod_Performance = BModEsts)
  
  return(results)
  
}


######Apply performance criteria#####
sim_ests <- sim_ests[lapply(sim_ests, length) > 0]
performance_results <- lapply(sim_ests, function(x) est_performance(x))
all_ests <- do.call(Map, c(f = rbind, performance_results))
save(all_ests, file = 'Estimator_Performance_Results_iNEXT.RData')


#-----Plot Performance Across All Repetitions-----
rm(list=ls()); gc()
load('Estimator_Performance_Results_iNEXT.RData')

#Plotting helpers
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

data_summary_0bound_pos <- function(x) {
  
  m <- mean(x)
  ymax <- m+sd(x)
  ymin <- m-sd(x)
  if(ymin <= 0){ymin <- 1e-5}
  return(c(y=m, ymin=ymin,ymax=ymax))
  
}

data_summary_0bound <- function(x) {
  
  m <- mean(x)
  ymax <- m+sd(x)
  if(ymax >= 0){ymax <- -1e-5}
  ymin <- m-sd(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
  
}

data_summary_01bound <- function(x) {
  
  m <- mean(x)
  ymin <- m-sd(x)
  if(ymin <= 0){ymin <- 1e-5}
  ymax <- m+sd(x)
  if(ymax >= 1){ymax <- 1-1e-5}
  return(c(y=m, ymin=ymin,ymax=ymax))
  
}


metrics <- c('Mean_Bias', 
             'Mean_Accuracy', 
             #'Mean_Coverage', 
             'Precision')

#####Plot functions#####
plotty <- function(x, type, estimators, metrics, scenarios, n_samples = c(3, 6, 12),
                   text_size = 16, labels = c('Obs', 'Chao', 'MSOM'), sub_plot_margins = c(7.5, 10, 7.5, 10), y_margins = c(22.5, 17),
                   x_margins = c(3,3), y_title = rep(T, length(scenarios)), x_labels = rep(T, length(scenarios)),
                   acc_breaks = waiver(), bias_breaks = waiver(), prec_breaks = waiver(), 
                   cov_breaks = seq(0, 1, 0.2), w95_breaks = waiver()){
  
  #Plot holders
  sub_plot_list <- vector('list', length(metrics)) 
  names(sub_plot_list) <- metrics
  x_plots <- vector('list', length(scenarios))
  names(x_plots) <- scenarios
  x_plots <- lapply(x_plots, function(x) x_plots[x] <- sub_plot_list)
  
  #Plot labels
  nsamples_labels <- c(`3` = '3 Samples',
                       `6` = '6 Samples',
                       `12` = '12 Samples')
  nsamples_labels <- nsamples_labels[names(nsamples_labels) %in% n_samples]
  
  x <- x[x$Estimator %in% estimators & x$Detection_Scenario %in% scenarios & x$Sampling_Events %in% n_samples, c('Estimator', 'Sampling_Events', 'Detection_Scenario', metrics)]
  
  s_counter <- 1
  
  for(s in scenarios){
    for(m in metrics){     
      
      #Isolate plot data
      plot_data <- x[x$Detection_Scenario == s, c('Estimator', 'Sampling_Events', m)]
      names(plot_data)[3] <- 'Metric'
      
      #Create base plot
      x_plots[[as.character(s)]][[m]] <- 
        ggplot(plot_data,
               aes(y = Metric,
                   x = Estimator)) +
        geom_violin(trim = T, size = 1, bw = 0.04, 
                    colour = NA,  
                    fill = alpha('#DB0F14', 0.2)) +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_line(size = 1.5),
              plot.margin = unit(sub_plot_margins, 'pt'), 
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              panel.border = element_rect(size = 1.5))
      
      if(y_title[s_counter]){
       
         x_plots[[as.character(s)]][[m]] <- 
          x_plots[[as.character(s)]][[m]] + 
            ylab(switch(m, 
                        'Mean_Accuracy' = 'Accuracy',
                        'Mean_Bias' = 'Bias',
                        'Precision' = 'Precision',
                        'Mean_Coverage' = 'Coverage',
                        'Mean_Width_95' = 'Uncertainty')) 
        
      } else {
        
        x_plots[[as.character(s)]][[m]] <- 
          x_plots[[as.character(s)]][[m]] + 
          ylab('') 
        
      }
      
      if(x_labels[s_counter]){
        
        x_plots[[as.character(s)]][[m]] <- 
          x_plots[[as.character(s)]][[m]] + 
          scale_x_discrete(labels = labels) 
        
      } else {
        
        x_plots[[as.character(s)]][[m]] <- 
          x_plots[[as.character(s)]][[m]] + 
          scale_x_discrete(labels = NULL) + 
          theme(axis.text.x = element_blank())
        
        
      }
      
      #Facet Plots and set appropriate axes for each metric
      if(m == 'Mean_Bias'){
        
        x_plots[[as.character(s)]][[m]] <-
          x_plots[[as.character(s)]][[m]] +
          scale_y_continuous(limits = c(floor_dec(quantile(x$Mean_Bias, 0.025, na.rm = T)) - 0.2,
                                        ceiling_dec(quantile(x$Mean_Bias, 0.975, na.rm = T)) + 0.05),
                             breaks = bias_breaks,
                             labels = scales::label_number(0.01)) +
          geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed', size = 1) +
          stat_summary(fun.data = mean_sd,
                       geom = "pointrange", color = "black", size = 0.5) +
          stat_summary(fun.data = mean_sd,
                       geom = "errorbar", color = "black", size = 0.5, width = 0.1) +
          theme(axis.title.y = element_text(face = 'bold', vjust = 1.5, size = text_size),
                axis.text.y = element_text(face = 'bold', size = text_size))
        
        if(length(n_samples) > 1){
          
          x_plots[[as.character(s)]][[m]] <-
            x_plots[[as.character(s)]][[m]] +
            facet_grid(~Sampling_Events, 
                       labeller = as_labeller(nsamples_labels)) +
            theme(strip.text = element_text(face = 'bold', 
                                            size = text_size),
                  strip.background = element_rect(size = 2))
          
        }
        
      } else if(m == 'Mean_Accuracy'){
        
        x_plots[[as.character(s)]][[m]] <-
          x_plots[[as.character(s)]][[m]] +
          scale_y_continuous(limits = c(floor_dec(quantile(x$Mean_Accuracy, 0.025, na.rm = T)) - 0.1, 0.001),
                             breaks = acc_breaks,
                             labels = scales::label_number(0.01)) +
          stat_summary(fun.data = data_summary_0bound,
                       geom = "pointrange", color = "black", size = 0.5) +
          stat_summary(fun.data = data_summary_0bound,
                       geom = "errorbar", color = "black", size = 0.5, width = 0.1) +
          theme(axis.title.y = element_text(face = 'bold', vjust = 1.5, size = text_size),
                axis.text.y = element_text(face = 'bold', size = text_size))
        
        if(type == 'Mod_Ov'){
          
          x_plots[[as.character(s)]][[m]] <-
              x_plots[[as.character(s)]][[m]] +
              facet_grid(~Sampling_Events, 
                         labeller = as_labeller(nsamples_labels)) +
              theme(strip.text = element_text(face = 'bold', 
                                              size = text_size),
                                 strip.background = element_rect(size = 2)) 
        } else if(length(n_samples) > 1){
          
          x_plots[[as.character(s)]][[m]] <-
            x_plots[[as.character(s)]][[m]] +
            facet_grid(~Sampling_Events, 
                       labeller = as_labeller(nsamples_labels)) +
            theme(strip.text = element_blank(),
                  strip.background = element_blank()) 
          
        }
        
      } else if(m == 'Mean_Coverage'){
        
        x_plots[[as.character(s)]][[m]] <-
          x_plots[[as.character(s)]][[m]] +
          scale_y_continuous(limits = c(-0.02, 1.02),
                             breaks = cov_breaks,
                             labels = paste0('  ', seq(0,1,.2))) +
          stat_summary(fun.data = data_summary_01bound,
                       geom = "pointrange", color = "black", size = 0.5) +
          stat_summary(fun.data = data_summary_01bound,
                       geom = "errorbar", color = "black", size = 0.5, width = 0.1) +
          theme(axis.title.y = element_text(face = 'bold', vjust = 1.5, size = text_size),
                axis.text.y = element_text(face = 'bold', size = text_size))
        
        if(length(n_samples) > 1 & type != 'Uncertainty'){
          
          x_plots[[as.character(s)]][[m]] <-
            x_plots[[as.character(s)]][[m]] +
            facet_grid(~Sampling_Events, 
                       labeller = as_labeller(nsamples_labels)) +
            theme(strip.text = element_blank(),
                  strip.background = element_blank()) 
          
        }
        
        if(length(n_samples) > 1 & type == 'Uncertainty'){
          
          x_plots[[as.character(s)]][[m]] <-
            x_plots[[as.character(s)]][[m]] +
            facet_grid(~Sampling_Events, 
                       labeller = as_labeller(nsamples_labels)) +
            theme(strip.text = element_text(face = 'bold', 
                                            size = text_size),
                  strip.background = element_rect(size = 2)) 
          
        }
        
        
      } else if(m == 'Precision'){
        
        x_plots[[as.character(s)]][[m]] <-
          x_plots[[as.character(s)]][[m]] +
          scale_y_continuous(limits = c(floor_dec(quantile(x$Precision, 0.025, na.rm = T)), 0.001),
                             breaks = prec_breaks,
                             labels = scales::label_number(0.01)) +
          theme(axis.text.x = element_text(face = 'bold', size = text_size, angle = 25)) +
          stat_summary(fun.data = data_summary_0bound,
                       geom = "pointrange", color = "black", size = 0.5) +
          stat_summary(fun.data = data_summary_0bound,
                       geom = "errorbar", color = "black", size = 0.5, width = 0.1) +
          theme(axis.title.y = element_text(face = 'bold', vjust = 1.5, size = text_size),
                axis.text.y = element_text(face = 'bold', size = text_size))
        
        if(length(n_samples) > 1){
          
          x_plots[[as.character(s)]][[m]] <-
            x_plots[[as.character(s)]][[m]] +
            facet_grid(~Sampling_Events, 
                       labeller = as_labeller(nsamples_labels)) +
            theme(strip.text = element_blank(),
                  strip.background = element_blank()) 
          
        }
        
      } else if(m == 'Mean_Width_95'){
        
        x_plots[[as.character(s)]][[m]] <-
          x_plots[[as.character(s)]][[m]] +
          scale_y_continuous(limits = c(-0.001, ceiling_dec(quantile(x$Mean_Width_95, 0.975, na.rm = T)) + 0.2),
                             breaks = w95_breaks,
                             labels = scales::label_number(0.01)) +
          theme(axis.text.x = element_text(face = 'bold', size = text_size, angle = 25)) +
          stat_summary(fun.data = data_summary_0bound_pos,
                       geom = "pointrange", color = "black", size = 0.5) +
          stat_summary(fun.data = data_summary_0bound_pos,
                       geom = "errorbar", color = "black", size = 0.5, width = 0.1) +
          
          theme(axis.title.y = element_text(face = 'bold', vjust = 1.5, size = text_size),
                axis.text.y = element_text(face = 'bold', size = text_size))
        
        if(length(n_samples) > 1){
          
          x_plots[[as.character(s)]][[m]] <-
            x_plots[[as.character(s)]][[m]] +
            facet_grid(~Sampling_Events, 
                       labeller = as_labeller(nsamples_labels)) +
            theme(strip.text = element_blank(),
                  strip.background = element_blank()) 
          
        }
        
      }
    }
    
    s_counter <- s_counter + 1
    
  }
  
  #Combine metric plots for each detection scenario
  combined_x_plots <- vector('list', length(scenarios))
  
  if(type == 'Full'){
    for(i in 1:length(scenarios)){
      if(!i %% 3){
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1.1,1,#1,
                                                       1.25), align = 'v')
        
      } else { 
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                            ncol = 1, nrow = length(metrics),
                                            heights = c(1.2,1,#1,
                                                        1), align = 'v')
        
      }
      
      combined_x_plots[[i]]$theme$plot.margin <- unit(c(35, 15, 35, 15), 'pt')
      
    }
  } 
  
  if(type == 'Mod_Ov'){
    for(i in 1:length(scenarios)){
      
      x_plots[[i]][["Mean_Accuracy"]]$theme$plot.margin <- unit(c(22.5, 3, 22.5, 3), 'pt')
      
      if(!i %% 3){
        
        x_plots[[i]][["Mean_Accuracy"]] <- 
          x_plots[[i]][["Mean_Accuracy"]] + 
          theme(axis.text.x = element_text(face = 'bold', size = text_size, angle = 25))
        
      } else {
        
        x_plots[[i]][["Mean_Accuracy"]]$theme$axis.text.x <- element_blank()
        
      }
      
      combined_x_plots[[i]] <- x_plots[[i]][['Mean_Accuracy']]
      
    }
  } 
  
  if(type == 'Subset'){
    for(i in 1:length(scenarios)){
      
      if(!i %% 3){
        
        x_plots[[i]][[length(x_plots[[i]])]] <- 
          x_plots[[i]][[length(x_plots[[i]])]] + 
          theme(axis.text.x = element_text(face = 'bold', size = text_size, angle = 25))
        
        for(n in 1:(length(x_plots[[i]]) - 1)){
          
          x_plots[[i]][[n]]$theme$axis.text.x <- element_blank()
        
          }
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1.175,1,1.15), align = 'v')
        
      } else {
        
        for(n in 1:length(x_plots[[i]])){
          
          x_plots[[i]][[n]]$theme$axis.text.x <- element_blank()
          
        }
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1.18,1,1), align = 'v')
        
      }
      
      
      combined_x_plots[[i]]$theme$plot.margin <- unit(c(7.5, 3, 7.5, 3), 'pt')
      
    }
  } 
  
  if(type == 'Main_Text'){
    for(i in 1:length(scenarios)){
      if(i == 3){
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1,1,1.3), align = 'v') 
        
      } else {
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1,1,1), align = 'v') 
        
      }
      
      combined_x_plots[[i]]$theme$plot.margin <- unit(c(y_margins[1], x_margins[1], y_margins[2], x_margins[2]), 'pt')
      
    }
  } 
  
  if(type == 'Uncertainty'){
    
    for(i in 1:length(scenarios)){
      if(!i %% 3){
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1.1,1.2), align = 'v')
        
      } else { 
        
        combined_x_plots[[i]] <- ggarrange(plotlist = x_plots[[i]],
                                           ncol = 1, nrow = length(metrics),
                                           heights = c(1.1,1), align = 'v')
        
      }
      
      combined_x_plots[[i]]$theme$plot.margin <- unit(c(22.5, 3, 22.5, 3), 'pt')
      
    }
    
  }
  
  return(combined_x_plots)
  
}


#####Site-Level Richness#####
siteR <- all_ests$Richness_Performance
siteR <- siteR %>% ungroup() %>% 
  select(Patch, Area, True_R, Estimator,
         Mean_Bias, Mean_Accuracy, Mean_Coverage, 
         Precision, Sampling_Events,
         bArea_Shape, Detection_Scenario, ID) %>% 
  distinct()
siteR$Precision[is.na(siteR$Precision)] <- 0

#Factorise estimators to help plot discretely 
estimators <- c('Observed', 'Chao', 'HMSOM')
siteR$Estimator <- factor(siteR$Estimator, levels = estimators)

#Plot
combined_SiteR_plots <- plotty(siteR, 
                               type = 'Full',
                               estimators = estimators,
                               metrics = metrics,
                               scenarios = 1:9,
                               text_size = 30,
                               x_labels = c(0,0,1,0,0,1,0,0,1),
                               acc_breaks = c(0, -0.35, -0.7, -1.05),
                               bias_breaks = c(0, -0.4, -0.8),
                               prec_breaks = c(0, -0.25, -0.5))

#Combine all plots and save as png
hjusts <- c(-.925, -.975, -.975, -.975, -1.05, -1.05, -1.05, -1.125, -1.125)

SiteR_full_plot <- fs::path('./Figures/Site_Richness/',  "Site_Level_Richness.png")
agg_png(SiteR_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(combined_SiteR_plots[[1]], combined_SiteR_plots[[4]], combined_SiteR_plots[[7]],
          combined_SiteR_plots[[2]], combined_SiteR_plots[[5]], combined_SiteR_plots[[8]],
          combined_SiteR_plots[[3]], combined_SiteR_plots[[6]], combined_SiteR_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SiteR_full_plot)

######Site-Level Beta#####
siteB <- all_ests$Beta_Performance
siteB <- siteB %>% 
  ungroup() %>% 
  select(Site_1, Site_2, Estimator, Samples, True_B,  
         Mean_Bias, Mean_Accuracy, Mean_Coverage, Precision, 
         ID, Sampling_Events, Detection_Scenario) %>% 
  distinct()
siteB$Precision[is.na(siteB$Precision)] <- 0

b_estimators <- c('Observed', 'Chao', 'HMSOM')
siteB$Estimator <- factor(siteB$Estimator, levels = b_estimators)

combined_SiteB_plots <- plotty(siteB, 
                               type = 'Full', 
                               estimators = b_estimators,
                               metrics = metrics,
                               scenarios = 1:9,
                               text_size = 30,
                               x_labels = c(0,0,1,0,0,1,0,0,1),
                               acc_breaks = c(0, -0.15, -0.3),
                               bias_breaks = c(0.2, 0, -0.2, -0.4),
                               prec_breaks = c(0, -0.2, -0.4, -0.6))

SiteB_full_plot <- fs::path('./Figures/Pairwise_Beta/',  "Site_Level_Beta.png")
agg_png(SiteB_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(combined_SiteB_plots[[1]], combined_SiteB_plots[[4]], combined_SiteB_plots[[7]],
          combined_SiteB_plots[[2]], combined_SiteB_plots[[5]], combined_SiteB_plots[[8]],
          combined_SiteB_plots[[3]], combined_SiteB_plots[[6]], combined_SiteB_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SiteB_full_plot)

#####SAR Slope#####
SARz <- all_ests$SAR_Performance
SARz <- SARz %>% 
  ungroup() %>%
  select(Estimator, z_Mean_Bias, z_Mean_Accuracy, z_Mean_Coverage, z_Precision,
         ID, Detection_Scenario, Sampling_Events) %>%
  distinct()
names(SARz)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')

SARz$Estimator <- factor(SARz$Estimator, levels = estimators)

combined_SARz_plots <- plotty(SARz,  
                              type = 'Full',
                              estimators = estimators,
                              metrics = metrics,
                              scenarios = 1:9,
                              text_size = 30,
                              x_labels = c(0,0,1,0,0,1,0,0,1),
                              acc_breaks = c(0, -0.15, -0.3),
                              bias_breaks = c(0, -0.2, -0.4),
                              prec_breaks = c(0, -0.5, -1))

SARz_full_plot <- fs::path('./Figures/SAR/',  "SAR_Slope.png")
agg_png(SARz_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(combined_SARz_plots[[1]], combined_SARz_plots[[4]], combined_SARz_plots[[7]],
          combined_SARz_plots[[2]], combined_SARz_plots[[5]], combined_SARz_plots[[8]],
          combined_SARz_plots[[3]], combined_SARz_plots[[6]], combined_SARz_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SARz_full_plot)

#####Beta Model Slope#####
Beta_z <- all_ests$BetaMod_Performance
Beta_z$b_Precision[is.na(Beta_z$b_Precision)] <- 0
Beta_z <- Beta_z %>% 
  ungroup() %>%
  select(Estimator, b_Mean_Bias, b_Mean_Accuracy, b_Mean_Coverage, b_Precision,
         ID, Detection_Scenario, Sampling_Events) %>%
  distinct()
names(Beta_z)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')
Beta_z <- Beta_z[!is.na(Beta_z$Mean_Bias),]

Beta_z$Estimator <- factor(Beta_z$Estimator, levels = b_estimators)

combined_Beta_z_plots <- plotty(Beta_z,  
                                type = 'Full',
                                estimators = b_estimators,
                                metrics = metrics,
                                scenarios = 1:9,
                                text_size = 30,
                                x_labels = c(0,0,1,0,0,1,0,0,1),
                                acc_breaks = c(0, -0.15, -0.3),
                                bias_breaks = c(0.2, 0, -0.2),
                                prec_breaks = c(0, -0.15, -0.3))

Betaz_full_plot <- fs::path('./Figures/Beta_Model/',  "Beta_Model_Slope.png")
agg_png(Betaz_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(combined_Beta_z_plots[[1]], combined_Beta_z_plots[[4]], combined_Beta_z_plots[[7]],
          combined_Beta_z_plots[[2]], combined_Beta_z_plots[[5]], combined_Beta_z_plots[[8]],
          combined_Beta_z_plots[[3]], combined_Beta_z_plots[[6]], combined_Beta_z_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(Betaz_full_plot)

#####SAR Intercept#####
SARc <- all_ests$SAR_Performance
SARc <- SARc %>% 
  ungroup() %>%
  select(Estimator, c_Mean_Bias, c_Mean_Accuracy, c_Mean_Coverage, c_Precision,
         ID, Detection_Scenario, Sampling_Events) %>%
  distinct()
names(SARc)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')

SARc$Estimator <- factor(SARc$Estimator, levels = estimators)

combined_SARc_plots <- plotty(SARc, 
                              type = 'Full', 
                              estimators = estimators,
                              metrics = metrics,
                              scenarios = 1:9,
                              text_size = 30,
                              x_labels = c(0,0,1,0,0,1,0,0,1),
                              acc_breaks = c(0, -0.75, -1.5),
                              bias_breaks = c(1.5, 0.75, 0, -0.75),
                              prec_breaks = c(0, -0.15, -0.3))

SARc_full_plot <- fs::path('./Figures/SAR/',  "SAR_Intercept.png")
agg_png(SARc_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(combined_SARc_plots[[1]], combined_SARc_plots[[4]], combined_SARc_plots[[7]],
          combined_SARc_plots[[2]], combined_SARc_plots[[5]], combined_SARc_plots[[8]],
          combined_SARc_plots[[3]], combined_SARc_plots[[6]], combined_SARc_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SARc_full_plot)

#####Beta Model Intercept#####
Betac <- all_ests$BetaMod_Performance
Betac <- Betac %>% 
  ungroup() %>%
  select(Estimator, a_Mean_Bias, a_Mean_Accuracy, a_Mean_Coverage, a_Precision,
         ID, Detection_Scenario, Sampling_Events) %>%
  distinct()
names(Betac)[2:5] <- c('Mean_Bias', 'Mean_Accuracy', 'Mean_Coverage', 'Precision')
Betac <- Betac[!is.na(Betac$Mean_Bias),]

Betac$Estimator <- factor(Betac$Estimator, levels = b_estimators)

combined_Beta_c_plots <- plotty(Betac,  
                                type = 'Full',
                                estimators = b_estimators,
                                metrics = metrics,
                                scenarios = 1:9,
                                text_size = 30,
                                x_labels = c(0,0,1,0,0,1,0,0,1),
                                acc_breaks = c(0, -1, -2),
                                bias_breaks = c(1, 0, -1, -2),
                                prec_breaks = c(0, -.15, -.3))

Betac_full_plot <- fs::path('./Figures/Beta_Model/',  "Beta_Intercept.png")
agg_png(Betac_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(combined_Beta_c_plots[[1]], combined_Beta_c_plots[[4]], combined_Beta_c_plots[[7]],
          combined_Beta_c_plots[[2]], combined_Beta_c_plots[[5]], combined_Beta_c_plots[[8]],
          combined_Beta_c_plots[[3]], combined_Beta_c_plots[[6]], combined_Beta_c_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(Betac_full_plot)

#####Combined intercept and slope accuracy (overall model fit)#####
#SARs
SARov <- all_ests$SAR_Performance
SARov <- SARov %>%
  ungroup() %>% 
  select(Estimator, Mean_Model_Accuracy, ID, Detection_Scenario, Sampling_Events) %>% 
  distinct()
names(SARov)[2] <- 'Mean_Accuracy'

SARov$Estimator <- factor(SARov$Estimator, levels = estimators)

SARov_plots <- plotty(SARov,
                      type = 'Mod_Ov',
                      estimators = estimators,
                      metrics = 'Mean_Accuracy',
                      scenarios = 1:9,
                      text_size = 30,
                      x_labels = c(0,0,1,0,0,1,0,0,1))

SARov_full_plot <- fs::path('./Figures/SAR/',  "SAR_Overall_Accuracy.png")
agg_png(SARov_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(SARov_plots[[1]], SARov_plots[[4]], SARov_plots[[7]],
          SARov_plots[[2]], SARov_plots[[5]], SARov_plots[[8]],
          SARov_plots[[3]], SARov_plots[[6]], SARov_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SARov_full_plot)

#Beta models
Bov <- all_ests$BetaMod_Performance
Bov <- Bov%>%
  ungroup() %>% 
  select(Estimator, Mean_Model_Accuracy, ID, Detection_Scenario, Sampling_Events) %>% 
  distinct()
names(Bov)[2] <- 'Mean_Accuracy'

Bov$Estimator <- factor(Bov$Estimator, levels = b_estimators)

Bov_plots <- plotty(Bov,
                    type = 'Mod_Ov',
                    estimators = b_estimators,
                    metrics = 'Mean_Accuracy',
                    scenarios = 1:9,
                    text_size = 30,
                    x_labels = c(0,0,1,0,0,1,0,0,1))

Bov_full_plot <- fs::path('./Figures/Beta_Model/',  "Beta_Overall_Accuracy.png")
agg_png(Bov_full_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(Bov_plots[[1]], Bov_plots[[4]], Bov_plots[[7]],
          Bov_plots[[2]], Bov_plots[[5]], Bov_plots[[8]],
          Bov_plots[[3]], Bov_plots[[6]], Bov_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')


invisible(dev.off())
knitr::include_graphics(Bov_full_plot)


#####Uncertainty#####

#Richness
REsts <- all_ests$Richness_Performance

REsts <- REsts %>% 
  ungroup() %>% 
  select(Estimator, Samples, Mean_Coverage, Mean_Width_95,
         ID, Sampling_Events, Detection_Scenario) %>% 
  distinct() 
REsts$Estimator <- factor(REsts$Estimator, levels = c('Chao', 'HMSOM'))

r_uncert_plots <- plotty(REsts,
                         type = 'Uncertainty',
                         estimators = c('Chao', 'HMSOM'),
                         labels = c('Chao', 'MSOM'),
                         metrics = c('Mean_Coverage', 'Mean_Width_95'),
                         scenarios = 1:9,
                         text_size = 30,
                         x_labels = c(0,0,1,0,0,1,0,0,1))

SiteR_uncert_plot <- fs::path('./Figures/Site_Richness',  "Richness_Uncertainty.png")
agg_png(SiteR_uncert_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(r_uncert_plots[[1]], r_uncert_plots[[4]], r_uncert_plots[[7]],
          r_uncert_plots[[2]], r_uncert_plots[[5]], r_uncert_plots[[8]],
          r_uncert_plots[[3]], r_uncert_plots[[6]], r_uncert_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SiteR_uncert_plot)

#Beta
BEsts <- all_ests$Beta_Performance

BEsts <- BEsts %>% 
  ungroup() %>% 
  select(Estimator, Samples, Mean_Coverage, Mean_Width95,
         ID, Sampling_Events, Detection_Scenario) %>% 
  distinct() 
BEsts$Estimator <- factor(BEsts$Estimator, levels = c('Chao', 'HMSOM'))
names(BEsts)[4] <- 'Mean_Width_95'

b_uncert_plots <- plotty(BEsts,
                         type = 'Uncertainty',
                         estimators = c('Chao', 'HMSOM'),
                         labels = c('Chao', 'MSOM'),
                         metrics = c('Mean_Coverage', 'Mean_Width_95'),
                         scenarios = 1:9,
                         text_size = 30,
                         x_labels = c(0,0,1,0,0,1,0,0,1))

SiteB_uncert_plot <- fs::path('./Figures/Pairwise_Beta/',  "Beta_Uncertainty.png")
agg_png(SiteB_uncert_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(b_uncert_plots[[1]], b_uncert_plots[[4]], b_uncert_plots[[7]],
          b_uncert_plots[[2]], b_uncert_plots[[5]], b_uncert_plots[[8]],
          b_uncert_plots[[3]], b_uncert_plots[[6]], b_uncert_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SiteB_uncert_plot)

#SARs
SAREsts <- all_ests$SAR_Performance

SAREsts <- SAREsts %>% 
  ungroup() %>% 
  select(Estimator, z_Mean_Coverage, z_Mean_Width95,
         c_Mean_Coverage, c_Mean_Width95, ID, 
         Sampling_Events, Detection_Scenario) %>% 
  distinct()
SAREsts$Estimator <- factor(SAREsts$Estimator, 
                            levels = c('Observed', 'Chao', 'HMSOM'))

SARz <- SAREsts %>% select(Estimator, z_Mean_Coverage, z_Mean_Width95,
                           Sampling_Events, Detection_Scenario)
SARc <- SAREsts %>% select(Estimator, c_Mean_Coverage, c_Mean_Width95,
                           Sampling_Events, Detection_Scenario)
names(SARz)[2:3] <- names(SARc)[2:3] <- c('Mean_Coverage', 'Mean_Width_95')

z_uncert_plots <- plotty(SARz,
                         type = 'Uncertainty',
                         estimators = c('Observed', 'Chao', 'HMSOM'),
                         labels = c('Obs', 'Chao', 'MSOM'),
                         metrics = c('Mean_Coverage', 'Mean_Width_95'),
                         scenarios = 1:9,
                         text_size = 30,
                         x_labels = c(0,0,1,0,0,1,0,0,1))
c_uncert_plots <- plotty(SARc,
                         type = 'Uncertainty',
                         estimators = c('Observed', 'Chao', 'HMSOM'),
                         labels = c('Obs', 'Chao', 'MSOM'),
                         metrics = c('Mean_Coverage', 'Mean_Width_95'),
                         scenarios = 1:9,
                         text_size = 30,
                         x_labels = c(0,0,1,0,0,1,0,0,1))

SARz_uncert_plot <- fs::path('./Figures/SAR/',  "SAR_Slope_Uncertainty.png")
agg_png(SARz_uncert_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(z_uncert_plots[[1]], z_uncert_plots[[4]], z_uncert_plots[[7]],
          z_uncert_plots[[2]], z_uncert_plots[[5]], z_uncert_plots[[8]],
          z_uncert_plots[[3]], z_uncert_plots[[6]], z_uncert_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SARz_uncert_plot)

SARc_uncert_plot <- fs::path('./Figures/SAR/',  "SAR_Intercept_Uncertainty.png")
agg_png(SARc_uncert_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(c_uncert_plots[[1]], c_uncert_plots[[4]], c_uncert_plots[[7]],
          c_uncert_plots[[2]], c_uncert_plots[[5]], c_uncert_plots[[8]],
          c_uncert_plots[[3]], c_uncert_plots[[6]], c_uncert_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(SARc_uncert_plot)


#Beta Mods
BModEsts <- all_ests$BetaMod_Performance


BModEsts <- BModEsts %>% 
  ungroup() %>% 
  select(Estimator, b_Mean_Coverage, b_Mean_Width95,
         a_Mean_Coverage, a_Mean_Width95, ID, 
         Sampling_Events, Detection_Scenario) %>% 
  distinct() 
BModEsts$Estimator <- factor(BModEsts$Estimator, 
                             levels = c('Observed', 'Chao', 'HMSOM'))

BModb <- BModEsts %>% select(Estimator, b_Mean_Coverage, b_Mean_Width95,
                             Sampling_Events, Detection_Scenario)
BModa <- BModEsts %>% select(Estimator, a_Mean_Coverage, a_Mean_Width95,
                             Sampling_Events, Detection_Scenario)
names(BModa)[2:3] <- names(BModb)[2:3] <- c('Mean_Coverage', 'Mean_Width_95')

b_uncert_plots <- plotty(BModb,
                         type = 'Uncertainty',
                         estimators = c('Observed', 'Chao', 'HMSOM'),
                         labels = c('Obs', 'Chao', 'MSOM'),
                         metrics = c('Mean_Coverage', 'Mean_Width_95'),
                         scenarios = 1:9,
                         text_size = 30,
                         x_labels = c(0,0,1,0,0,1,0,0,1))
a_uncert_plots <- plotty(BModa,
                         type = 'Uncertainty',
                         estimators = c('Observed', 'Chao', 'HMSOM'),
                         labels = c('Obs', 'Chao', 'MSOM'),
                         metrics = c('Mean_Coverage', 'Mean_Width_95'),
                         scenarios = 1:9,
                         text_size = 30,
                         x_labels = c(0,0,1,0,0,1,0,0,1))

BModb_uncert_plot <- fs::path('./Figures/Beta_Model/',  "BetaMod_Slope_Uncertainty.png")
agg_png(BModb_uncert_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(b_uncert_plots[[1]], b_uncert_plots[[4]], b_uncert_plots[[7]],
          b_uncert_plots[[2]], b_uncert_plots[[5]], b_uncert_plots[[8]],
          b_uncert_plots[[3]], b_uncert_plots[[6]], b_uncert_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(BModb_uncert_plot)

BModa_uncert_plot <- fs::path('./Figures/Beta_Model/',  "BetaMod_Intercept_Uncertainty.png")
agg_png(BModa_uncert_plot, width = 1250, height = 800, units = "mm", res = 500) 

ggarrange(a_uncert_plots[[1]], a_uncert_plots[[4]], a_uncert_plots[[7]],
          a_uncert_plots[[2]], a_uncert_plots[[5]], a_uncert_plots[[8]],
          a_uncert_plots[[3]], a_uncert_plots[[6]], a_uncert_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 32, face = 'bold.italic'),
          hjust = hjusts, vjust = 1.1, heights = c(1,1,1.115),
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(BModa_uncert_plot)


#-----Subset plots-----
subset_siteR <- plotty(siteR, 
                       type = 'Subset',
                       estimators = estimators,
                       metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                       scenarios = c(4,5,6),
                       text_size = 16)
SiteR_sub_plot <- ggarrange(subset_siteR[[1]], subset_siteR[[2]], subset_siteR[[3]],
                             # labels = c('\u03bc = 0.005, \u03c3 = 1.0',
                             #            '\u03bc = 0.05, \u03c3 = 1.0',
                             #            '\u03bc = 0.2, \u03c3 = 1.0'),
                             # font.label = list(size = 26, face = 'bold.italic'),
                             # hjust = -1, vjust = 2,
                             ncol = 1, nrow = 3, heights = c(1,1,1.1))

subset_SARz <- plotty(SARz,  
                      type = 'Subset',
                      estimators = estimators,
                      metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                      scenarios = c(4,5,6),
                      text_size = 16)
SARz_sub_plot <- ggarrange(subset_SARz[[1]], subset_SARz[[2]], subset_SARz[[3]],
                            # labels = c('\u03bc = 0.005, \u03c3 = 1.0',
                            #            '\u03bc = 0.05, \u03c3 = 1.0',
                            #            '\u03bc = 0.2, \u03c3 = 1.0'),
                            # font.label = list(size = 26, face = 'bold.italic'),
                            # hjust = -1, vjust = 2,
                            ncol = 1, nrow = 3, heights = c(1,1,1.1))

subset_SARc <- plotty(SARc,  
                      type = 'Subset',
                      estimators = estimators,
                      metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                      scenarios = c(4,5,6),
                      text_size = 16)
SARc_sub_plot <- ggarrange(subset_SARc[[1]], subset_SARc[[2]], subset_SARc[[3]],
                            # labels = c('\u03bc = 0.005, \u03c3 = 1.0',
                            #            '\u03bc = 0.05, \u03c3 = 1.0',
                            #            '\u03bc = 0.2, \u03c3 = 1.0'),
                            # font.label = list(size = 26, face = 'bold.italic'),
                            # hjust = -1, vjust = 2,
                            ncol = 1, nrow = 3, heights = c(1,1,1.1))

rich_sub_plot <- fs::path('./Figures/Subsets/',  "Richness_subset.png")
agg_png(rich_sub_plot, width = 550, height = 500, units = "mm", res = 500) 

ggarrange(SiteR_sub_plot, SARz_sub_plot, SARc_sub_plot,
          ncol = 3, nrow = 1)

invisible(dev.off())
knitr::include_graphics(rich_sub_plot)


subset_siteB <- plotty(siteB,  
                       type = 'Subset',
                       estimators = estimators,
                       metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                       scenarios = c(4,5,6),
                       text_size = 16)
SiteB_sub_plot <- ggarrange(subset_siteB[[1]], subset_siteB[[2]], subset_siteB[[3]],
                             # labels = c('\u03bc = 0.005, \u03c3 = 1.0',
                             #            '\u03bc = 0.05, \u03c3 = 1.0',
                             #            '\u03bc = 0.2, \u03c3 = 1.0'),
                             # font.label = list(size = 26, face = 'bold.italic'),
                             # hjust = -1, vjust = 2,
                             ncol = 1, nrow = 3, heights = c(1,1,1.1))

subset_Beta_z <- plotty(Beta_z,  
                        type = 'Subset',
                        estimators = estimators,
                        metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                        scenarios = c(4,5,6),
                        text_size = 16)
Beta_z_sub_plot <- ggarrange(subset_Beta_z[[1]], subset_Beta_z[[2]], subset_Beta_z[[3]],
                              # labels = c('\u03bc = 0.005, \u03c3 = 1.0',
                              #            '\u03bc = 0.05, \u03c3 = 1.0',
                              #            '\u03bc = 0.2, \u03c3 = 1.0'),
                              # font.label = list(size = 26, face = 'bold.italic'),
                              # hjust = -1, vjust = 2,
                              ncol = 1, nrow = 3, heights = c(1,1,1.1))

subset_Betac <- plotty(Betac,  
                       type = 'Subset',
                       estimators = estimators,
                       metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                       scenarios = c(4,5,6),
                       text_size = 16)
Betac_sub_plot <- ggarrange(subset_Betac[[1]], subset_Betac[[2]], subset_Betac[[3]],
                             # labels = c('\u03bc = 0.005, \u03c3 = 1.0',
                             #            '\u03bc = 0.05, \u03c3 = 1.0',
                             #            '\u03bc = 0.2, \u03c3 = 1.0'),
                             # font.label = list(size = 26, face = 'bold.italic'),
                             # hjust = -1, vjust = 2,
                             ncol = 1, nrow = 3, heights = c(1,1,1.1))

b_sub_plot <- fs::path('./Figures/Subsets/',  "Beta_subset.png")
agg_png(b_sub_plot, width = 550, height = 500, units = "mm", res = 500) 

ggarrange(SiteB_sub_plot, Beta_z_sub_plot, Betac_sub_plot,
          ncol = 3, nrow = 1)

invisible(dev.off())
knitr::include_graphics(b_sub_plot)


#-----  Main Text Plots -----
#Richness
main_siteR <- plotty(siteR, 
                     type = 'Main_Text',
                     estimators = estimators,
                     metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                     scenarios = c(4,5,6),
                     n_samples = 6,
                     text_size = 18,
                     y_margins = c(10,10),
                     x_margins = c(1,1),
                     sub_plot_margins = c(5, 12.5, 5, 12.5),
                     y_title = c(T,T,T), x_labels = c(F,F,T),
                     bias_breaks = c(0, -0.4, -0.8),
                     acc_breaks = c(0, -0.25, -0.5, -0.75),
                     prec_breaks = c(0, -0.2, -0.4, -0.6))

main_SARz <- plotty(SARz,  
                    type = 'Main_Text',
                      estimators = estimators,
                      metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                      scenarios = c(4,5,6),
                      n_samples = 6,
                      text_size = 18,
                      y_margins = c(10,10),
                    x_margins = c(1,1),
                    sub_plot_margins = c(5, 12.5, 5, 12.5),
                    y_title = c(F,F,F), x_labels = c(F,F,T),
                    bias_breaks = c(0, -0.2, -0.4, -0.6),
                    acc_breaks = c(0, -0.15, -0.3),
                    prec_breaks = c(0, -0.5, -1))

main_SARc <- plotty(SARc,  
                    type = 'Main_Text',
                    estimators = estimators,
                    metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                    scenarios = c(4,5,6),
                    n_samples = 6,
                    text_size = 18,
                    y_margins = c(10,10),
                    x_margins = c(1,1),
                    sub_plot_margins = c(5, 12.5, 5, 12.5),
                    y_title = c(F,F,F), x_labels = c(F,F,T),
                    bias_breaks = c(1.5, 1, 0.5, 0),
                    acc_breaks = c(0, -0.4, -0.8, -1.2),
                    prec_breaks = c(0, -0.15, -0.3))


main_siteR_plot <- ggarrange(plotlist = main_siteR, 
                             ncol = 1, nrow = 3, align = 'hv', heights = c(1,1,1.15))
main_SARz_plot <- ggarrange(plotlist = main_SARz, 
                             ncol = 1, nrow = 3, align = 'hv', heights = c(1,1,1.15))
main_SARc_plot <- ggarrange(plotlist = main_SARc, 
                             ncol = 1, nrow = 3, align = 'hv', heights = c(1,1,1.15))

rich_main_plot <- fs::path('./Figures/Main/',  "Richness_Main.png")
agg_png(rich_main_plot, width = 350, height = 475, units = "mm", res = 500) 

ggarrange(main_siteR_plot, main_SARz_plot, main_SARc_plot,
          ncol = 3, nrow = 1, align = 'hv')

invisible(dev.off())
knitr::include_graphics(rich_main_plot)


#Beta
main_siteB <- plotty(siteB,  
                     type = 'Main_Text',
                     estimators = estimators,
                     metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                     scenarios = c(4,5,6),
                     n_samples = 6,
                     text_size = 18,
                     y_margins = c(10,10),
                     x_margins = c(1,1),
                     sub_plot_margins = c(5, 12.5, 5, 12.5),
                     y_title = c(T,T,T), x_labels = c(F,F,T),
                     bias_breaks = c(0.25, 0, -0.25, -0.5),
                     acc_breaks = c(0, -0.15, -0.3),
                     prec_breaks = c(0, -0.2, -0.4))

main_Beta_z <- plotty(Beta_z,  
                      type = 'Main_Text',
                      estimators = estimators,
                      metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                      scenarios = c(4,5,6),
                      n_samples = 6,
                      text_size = 18,
                      y_margins = c(10,10),
                      x_margins = c(1,1),
                      sub_plot_margins = c(5, 12.5, 5, 12.5),
                      y_title = c(F,F,F), x_labels = c(F,F,T),
                      bias_breaks = c(0.2, 0, -0.2),
                      acc_breaks = c(0, -0.15, -0.3),
                      prec_breaks = c(0, -0.2, -0.4)
                      )

main_Beta_c <- plotty(Betac,  
                      type = 'Main_Text',
                      estimators = estimators,
                      metrics = c('Mean_Bias', 'Mean_Accuracy', 'Precision'),
                      scenarios = c(4,5,6),
                      n_samples = 6,
                      text_size = 18,
                      y_margins = c(10,10),
                      x_margins = c(1,1),
                      sub_plot_margins = c(5, 12.5, 5, 12.5),
                      y_title = c(F,F,F), x_labels = c(F,F,T),
                      bias_breaks = c(1, 0, -1, -2),
                      acc_breaks = c(0, -1, -2),
                      prec_breaks = c(0, -0.15, -0.3)
                      )


main_siteB_plot <- ggarrange(plotlist = main_siteB, 
                             ncol = 1, nrow = 3, align = 'hv', heights = c(1,1,1.15))
main_Beta_z_plot <- ggarrange(plotlist = main_Beta_z, 
                              ncol = 1, nrow = 3, align = 'hv', heights = c(1,1,1.15))
main_Beta_c_plot <- ggarrange(plotlist = main_Beta_c, 
                              ncol = 1, nrow = 3, align = 'hv', heights = c(1,1,1.15))

beta_main_plot <- fs::path('./Figures/Main/',  "Beta_Main.png")
agg_png(beta_main_plot, width = 350, height = 475, units = "mm", res = 500) 

ggarrange(main_siteB_plot, main_Beta_z_plot, main_Beta_c_plot,
          ncol = 3, nrow = 1, align = 'hv')

invisible(dev.off())
knitr::include_graphics(beta_main_plot)


#-----Supplementary analyses and plots-----

#####Example SARs#####
REsts <- all_ests$Richness_Performance

REsts <- REsts %>% 
  ungroup() %>% 
  select(ID, Estimator, Patch, Area, Estimate, True_R, Samples, Repeat,
         Detection_Scenario, bArea_Shape, Sampling_Events) %>% 
  distinct() %>% 
  filter(Repeat == 1) %>% 
  filter(Sampling_Events == 6 & bArea_Shape == 8)

ids <- vector(length = 9)
for(i in 1:9){
  ds_sub <- REsts[REsts$Detection_Scenario == i,]
  ids[i] <- ds_sub$ID[sample(1:length(ds_sub$ID), 1)]
}

REsts <- all_ests$Richness_Performance
REsts <- REsts[REsts$ID %in% ids & REsts$Repeat == 1,]
REsts <- REsts %>% 
  filter(Estimator %in% c('Observed', 'Chao', 'HMSOM')) %>% 
  filter(Samples == 'Pooled' | is.na(Samples))
sars_list <- split(REsts, f = REsts$Detection_Scenario)
sar_plots <- vector('list', 9)

for(i in 1:9){
  
  true <- unique(sars_list[[i]][, c('Area', 'True_R')]) %>% na.omit()  
  names(true) <- c('Area', 'Estimate')
  true$Estimator <- 'True'
  true <- true[,c(1,3,2)]
  ests <- sars_list[[i]][,c('Area', 'Estimator', 'Estimate')]
  ests <- rbind(ests, true)
  ests$Estimator[ests$Estimator == 'HMSOM'] <- 'MSOM'
  ests$Estimator <- factor(ests$Estimator, levels = c('True', 'Observed', 'Chao', 'MSOM'))
  sar_plots[[i]] <- 
    ggplot(ests, aes(x = log(Area), y = log(Estimate), colour = Estimator, group = Estimator)) +
    geom_point(size = 3) +
    geom_smooth(aes(group = Estimator), formula = y ~ x, method = lm, se = F, size = 1.5) +
    scale_color_viridis_d() +
    scale_y_continuous(limits = log(c(min(REsts$Estimate), 200)), breaks = log(c(25, 50, 100, 200)), labels = c(25, 50, 100, 200)) +
    scale_x_continuous(limits = log(c(15, 20100)), breaks = log(c(100,1000,10000)), labels = c(100,1000,10000)) +
    labs(x = ''#bquote('\nPatch Area (Ha)')
         , y = ''#'Species Richness#\n'
         ) +
    theme_bw() +
    theme(axis.text = element_text(face = 'bold', colour = 'black', size = 16),
          #axis.title = element_text(face = 'bold', colour = 'black', size = 12),
          panel.grid = element_blank(),
          axis.line = element_line(size = 1.1),
          axis.ticks = element_line(size = 1.1),
          panel.border = element_blank(),
          plot.background = element_blank(),
          legend.position = 'none',
          plot.margin = unit(c(15, 3, 15, 3), 'pt'))
  
}

sar_plot <- fs::path('./Figures/Summaries/',  "Example_SARs.png")
agg_png(sar_plot, width = 500, height = 500, units = "mm", res = 500) 

ggarrange(sar_plots[[1]], sar_plots[[4]], sar_plots[[7]],
          sar_plots[[2]], sar_plots[[5]], sar_plots[[8]],
          sar_plots[[3]], sar_plots[[6]], sar_plots[[9]],
          labels = c('Scenario 1: \u03bc = 0.005, \u03c3 = 0.25',
                     'Scenario 2: \u03bc = 0.005, \u03c3 = 1.0',
                     'Scenario 3: \u03bc = 0.005, \u03c3 = 3.0',
                     'Scenario 4: \u03bc = 0.05, \u03c3 = 0.25',
                     'Scenario 5: \u03bc = 0.05, \u03c3 = 1.0',
                     'Scenario 6: \u03bc = 0.05, \u03c3 = 3.0',
                     'Scenario 7: \u03bc = 0.2, \u03c3 = 0.25',
                     'Scenario 8: \u03bc = 0.2, \u03c3 = 1.0',
                     'Scenario 9: \u03bc = 0.2, \u03c3 = 3.0'), 
          font.label = list(size = 22, face = 'bold.italic'),
          hjust = c(-.375, -.375, -.375, -.4, -.4, -.4, -.44, -.44, -.44),
          vjust = 1, legend = 'none',
          ncol = 3, nrow = 3, align = 'hv')

invisible(dev.off())
knitr::include_graphics(sar_plot)


leg <- ggplot(ests, aes(x = log(Area), y = log(Estimate), colour = Estimator, group = Estimator)) +
  geom_point(size = 3) +
  geom_smooth(aes(group = Estimator), formula = y ~ x, method = lm, se = F, size = 1.5) +
  scale_color_viridis_d(name = '') +
  scale_y_continuous(limits = log(c(min(REsts$Estimate), 200)), breaks = log(c(25, 50, 100, 200)), labels = c(25, 50, 100, 200)) +
  scale_x_continuous(limits = log(c(15, 20100)), breaks = log(c(100,1000,10000)), labels = c(100,1000,10000)) +
  labs(x = ''#bquote('\nPatch Area (Ha)')
       , y = ''#'Species Richness#\n'
  ) +
  theme_bw() +
  theme(axis.text = element_text(face = 'bold', colour = 'black', size = 16),
        #axis.title = element_text(face = 'bold', colour = 'black', size = 12),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.1),
        axis.ticks = element_line(size = 1.1),
        panel.border = element_blank(),
        plot.background = element_blank(),
        legend.position = 'top',
        legend.text= element_text(face = 'bold', size = 12),
        plot.margin = unit(c(15, 3, 15, 3), 'pt'))

leg_plot <- fs::path('./Figures/Summaries/',  "Legend.png")
agg_png(leg_plot, width = 250, height = 250, units = "mm", res = 500) 

leg

invisible(dev.off())
knitr::include_graphics(leg_plot)

#####Effect of area/difference in area on estimate bias#####

#Richness
REsts <- all_ests$Richness_Performance
REsts <- REsts %>% 
  ungroup() %>% 
  select(ID, Estimator, Patch, Area, Mean_Bias, 
         mu_ip, sd_ip, bArea_Shape, Sampling_Events) %>% 
  distinct()

obs_mod <- lm(Mean_Bias ~ mu_ip + sd_ip + Sampling_Events + log10(Area), 
                      data = REsts[REsts$Estimator == 'Observed',])
chao_mod <- lm(Mean_Bias ~ mu_ip + sd_ip + Sampling_Events + log10(Area), 
                       data = REsts[REsts$Estimator == 'Chao',])
msom_mod <- lm(Mean_Bias ~ mu_ip + sd_ip + Sampling_Events + log10(Area),
                       data = REsts[REsts$Estimator == 'HMSOM',])

summary(obs_mod)
summary(chao_mod)
summary(msom_mod)

obs_df <- ggpredict(obs_mod, terms = 'Area [25:20000]')
obs_df$Method <- 'Observed'
chao_df <- ggpredict(chao_mod, terms = 'Area [25:20000]')
chao_df$Method <- 'Chao'
msom_df <- ggpredict(msom_mod, terms = 'Area [25:20000]')
msom_df$Method <- 'MSOM'

r_bias_df <- rbind(obs_df, chao_df, msom_df)
r_bias_df$Method <- factor(r_bias_df$Method, levels = c('Observed', 'Chao', 'MSOM'))

r_bias_plot <- 
  ggplot(r_bias_df, aes(log10(x), predicted, 
                        colour = Method
                        )) +
      #geom_ribbon(aes(x = log10(x), ymin = conf.low, ymax = conf.high), alpha = 0) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed', size = 0.5) +
      labs(x = '\nRelative Patch Area', y = 'Bias\n') +
      scale_x_continuous(limits = c(log10(25), log10(21000)), 
                         labels = c(100,1000,10000), 
                         breaks = log10(c(100,1000,10000))) +
      scale_color_viridis_d() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(face = 'bold', size = 10, colour = 'black'),
            axis.title = element_text(face = 'bold', size = 11, colour = 'black'), 
            panel.border = element_blank(),
            axis.line = element_line(size = 1),
            axis.ticks = element_line(size = 1),
            legend.title = element_text(face = 'bold', size = 11, colour = 'black'), 
            legend.text = element_text(size = 10, colour = 'black'))  

r_bias_file <- fs::path('./Figures/Summaries/',  "Richness_Bias_v_Area.png")
agg_png(r_bias_file, width = 200, height = 150, units = "mm", res = 500) 

r_bias_plot

invisible(dev.off())
knitr::include_graphics(r_bias_file)

#Sorenson
areas <- REsts %>% select(ID, Patch, Area) %>% distinct()

BEsts <- all_ests$Beta_Performance
BEsts <- BEsts %>% 
  ungroup() %>% 
  select(ID, Estimator, Site_1, Site_2, Mean_Bias, 
         mu_ip, sd_ip, bArea_Shape, Sampling_Events) %>% 
  distinct()

uniq_pairs <- BEsts %>% select(ID, Site_1, Site_2) %>% distinct()
ids <- unique(uniq_pairs$ID)
pairs <- sum(uniq_pairs$ID == 1)
diffs <- vector(length = nrow(uniq_pairs))

for(i in 1:length(ids)){
  areas_ID <- areas[areas$ID == ids[i],]
  pairs_ID <- uniq_pairs[uniq_pairs$ID == ids[i],]
  for(n in 1:pairs){
    Area_1 <- areas_ID$Area[areas_ID$Patch == pairs_ID$Site_1[n]]
    Area_2 <- areas_ID$Area[areas_ID$Patch == pairs_ID$Site_2[n]]
    diffs[((i-1)*pairs)+n] <- abs(Area_1 - Area_2)
  }
}

uniq_pairs$Area_Diff <- diffs
BEsts <- merge(BEsts, uniq_pairs)

b_obs_mod <- lm(Mean_Bias ~ mu_ip + bArea_Shape + sd_ip + Sampling_Events + log10(Area_Diff + 1e-5), 
                   data = BEsts[BEsts$Estimator == 'Observed',])
b_chao_mod <- lm(Mean_Bias ~ mu_ip + bArea_Shape + sd_ip + Sampling_Events + log10(Area_Diff + 1e-5), 
                   data = BEsts[BEsts$Estimator == 'Chao',])
b_msom_mod <- lm(Mean_Bias ~ mu_ip + bArea_Shape + sd_ip + Sampling_Events + log10(Area_Diff + 1e-5), 
                   data = BEsts[BEsts$Estimator == 'HMSOM',])

summary(b_obs_mod)
summary(b_chao_mod)
summary(b_msom_mod)

b_obs_df <- ggpredict(b_obs_mod, terms = 'Area_Diff [1:19975]')
b_obs_df$Method <- 'Observed'
b_chao_df <- ggpredict(b_chao_mod, terms = 'Area_Diff [1:19975]')
b_chao_df$Method <- 'Chao'
b_msom_df <- ggpredict(b_msom_mod, terms = 'Area_Diff [1:19975]')
b_msom_df$Method <- 'MSOM'

b_bias_df <- rbind(b_obs_df, b_chao_df, b_msom_df)
b_bias_df$Method <- factor(b_bias_df$Method, levels = c('Observed', 'Chao', 'MSOM'))

b_bias_plot <- 
  ggplot(b_bias_df, aes(log10(x), predicted, 
                        colour = Method
                        )) +
      #geom_ribbon(aes(x = log10(x), ymin = conf.low, ymax = conf.high), alpha = 0) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed', size = 0.5) +
      labs(x = '\nRelative Difference in Patch Area', y = 'Bias\n') +
      scale_x_continuous(limits = c(log10(1), log10(20000)), 
                         labels = c(1,10,100,1000,10000), 
                         breaks = log10(c(1,10,100,1000,10000))) +
      scale_color_viridis_d() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(face = 'bold', size = 10, colour = 'black'),
            axis.title = element_text(face = 'bold', size = 11, colour = 'black'), 
            panel.border = element_blank(),
            axis.line = element_line(size = 1),
            axis.ticks = element_line(size = 1),
            legend.title = element_text(face = 'bold', size = 11, colour = 'black'), 
            legend.text = element_text(size = 10, colour = 'black'))  

b_bias_file <- fs::path('./Figures/Summaries/',  "Beta_Bias_v_Area.png")
agg_png(b_bias_file, width = 200, height = 150, units = "mm", res = 500) 

b_bias_plot

invisible(dev.off())
knitr::include_graphics(b_bias_file)