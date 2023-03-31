rm(list=ls()); gc()

library(pacman)
p_load(dplyr, mobsim, data.table, tidyr, doSNOW, parallel, sf, sp, betafunctions)

setwd()

#-----Simulation Parameters-----

n_patches <- 25  #Number of patches in each landscape
max_patch <- 20000  #Maximum patch size in Ha
min_patch <- 25 #Defined smallest patch size to give good range

patch_visits <- c(2, 3, 5)  #Number of visits to each patch

mu_ip <- c(0.025, 0.1, 0.2) #Mean individual level detection probability
sd_ip <- c(0.25, 2)  #SD of individual-level detection probabilities
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
params <- params[rep(seq_len(nrow(params)), 30), ] #Replicate param combos 30 times
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
  #Generate responses with four-parameter beta distribution
  spp_pot$alpha <- seq(0.1, bArea_Shape, length.out = nSpp) #Species level shape parameters, 
                                                            #increasing chance of severe area response with increasing rarity
                                                            #ensures realistic, viable patch-level populations of rarest species where they occur 
                                                            #i.e. if a species has a landscape level population of 500 individuals, they are likely to be concentrated in a few patches
  spp_pot$bArea <- round(mapply(FUN = rBeta.4P,                              
                                n = 1,
                                u = 5, 
                                l = 0, 
                                alpha = spp_pot$alpha,
                                beta = 5), 3)  #Generate species specific area responses through random draws from four-parameter beta distributions
  
  #Select 1/8th of species to have negative area responses
  nNegSpp <- round(nSpp/8)
  negSpp <- sample(1:nSpp, nNegSpp) #Randomly select 1/8th of species
  spp_pot$bArea[negSpp] <- -spp_pot$bArea[negSpp] #Assign generated negative response
  
  #Generate probabilities of assignment to each patch for each sp 
  assignment_probs <- matrix(NA, nrow = nSpp, ncol = n_patches)
  
  for(i in 1:nSpp){
    for(x in 1:n_patches){
      #Probability of assignment to patch = Exponent of area response * log(Patch Area)
      assignment_probs[i,x] <- exp(spp_pot$bArea[i] * log(patches[x]))
    }
    #Scale to sum to 1
    assignment_probs[i,] <- (assignment_probs[i,] / sum(assignment_probs[i,])) 
  }
  
  #Create empty Site x Species abundance matrix
  true_abun <- matrix(0, nrow = n_patches, ncol = nSpp)
  
  locs <- 1:n_patches #List of fragment IDs (used to remove fragments from pool once capacity reached)
  done <- FALSE
  
  while(sum(spp_pot$Count) > 0){ #While individuals in the species pool remain to be assigned...
    for(sp in 1:nSpp){ #Repeatedly assign 100 individuals from each species until all done
      if(spp_pot$Count[sp] > 0){ #If species individuals remain to be assigned proceed
        vacancies <- 0 #Set vacancies to 0 to enter next loop
        while(vacancies == 0){ #If there are 0 vacancies in patch selected in next section, retry assignment
          patch <- sample(locs, 1, prob = assignment_probs[sp, locs]) #Randomly select patch according to generated assignment probabilities
          vacancies <- patch_pops$Population[patch] #Find no. of individual vacancies yet to be filled in the patch
          if(vacancies > 0){ #If vacancies available in patch...
            true_abun[patch, sp] <- (true_abun[patch, sp] + 100) #Assign 100 individuals to the patch
            spp_pot$Count[sp] <- (spp_pot$Count[sp] - 100) #Decrease species count by 100
            patch_pops$Population[patch] <- (patch_pops$Population[patch] - 100) #Decrease remaining vacancies in patch by 100
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
    if(sum(spp_pot$Count) %% 10000 == 0){cat('\n', sum(spp_pot$Count))} #Counter in console
    if(done){break}
  }
  
  #Determine true richness at each site
  true_pa <- ifelse(true_abun > 0, 1, 0) #Convert to Presence absence matrix
  True_R <- rowSums(true_pa, na.rm = T)  #Summate spp richness
  patchR <- data.frame(Patch = 1:n_patches, Area = patches, True_R = True_R)
  #plot(patchR$Area, patchR$True_R, ylim = c(0, nSpp))
  #summary(lm(log(True_R) ~ log(Area), data = patchR))
  
  #Assess whether species density increases/decreases roughly linearly with area
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
  patch_proportion <- as.numeric(((search_area/patches) * n_transects)) #calculate proportion of each patch covered by transects
  
  #Species-specific, individual-level detection probabilities
  ip <- plogis(rnorm(nSpp, qlogis(mu_ip), sd_ip))  #(mu & sd on logit scale; back to probability)
  
  observed_data <- vector('list', 10)
  
  #Simulate 10 repeat sampling processes
  for(i in 1:10){
    
    set.seed(i)
    
    #Number of individuals within sampled area
    sampleable <- matrix(NA, nrow = n_patches, ncol = nSpp)
    for(j in 1:n_patches){
      for(x in 1:nSpp){
        sampleable[j,x] <- rbinom(1, true_abun[j,x], patch_proportion[j]) 
        #sampleable[j,x] <- round(true_abun[j,x] * patch_proportion[j])
      }
    }
    
    #Calculate sampleble richness
    #Sampleable <- rowSums(sampleable > 0)
    #plot(patchR$Area, Sampleable,  ylim = c(0, nSpp))
    #lm(log(Sampleable) ~ log(patchR$Area))
    
    #Sampling visits
    Y <- array(NA, dim = c(n_patches, nSpp, patch_visits))
    for(j in 1:n_patches){
      for(x in 1:nSpp){
        Y[j,x, ] <- rbinom(n = patch_visits, 
                           size = sampleable[j, x], 
                           prob = ip[x]) #Simulate imperfect detection using species individual-level detection probabilities
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
                                               'betafunctions', 'dqrng'),
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

save(landscape_simulations, file = 'Landscape_Community_Sims.RData')

######Parse simulation repeats into individual sampling processes for running on HPC cluster#####
landscape_sims_parsed <- vector('list', (nrow(params)*10))
x <- 1

for(i in 1:length(landscape_simulations)){
  for(n in seq(1, 10, 1)){
    
    landscape_sims_parsed[[x]] <- list(id = landscape_simulations[[i]][[1]],
                                       Parameters = landscape_simulations[[i]][[2]],
                                       True_Data = landscape_simulations[[i]][[3]],
                                       Observed_Data = list(
                                       landscape_simulations[[i]][[4]][[n]]))
    x <- x+1
  }
}

save(landscape_sims_parsed, file = 'Landscape_Sims_Parsed.RData')
