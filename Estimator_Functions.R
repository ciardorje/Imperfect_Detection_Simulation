rm(list=ls());gc()

library(pacman)
p_load(nimble, SpadeR, betapart, reshape2)

load('Landscape_Sims_Parsed.RData')

#####Estimation Functions#####

#Formula for sorensen index
sorensen <- function(x, y){ ((2 * sum(x * y)) / (sum(x) + sum(y))) }

#Probability complement function - calculates the probability that a species will be occurr at least once across n transect groups
#as a complement of the probability the species will not occur in any of the n transect groups
p_comp <- nimbleFunction(
  
  run = function(p = double(0), n = double(0)){
    
    pc <- (1 - p)
    pcn <- pc^n
    out <- (1 - pcn)
    
    return(out)
    returnType(double(0))
    
  })

#Hierarchical Multi-Species Occupancy Model structure
HMSOM <- nimbleCode({
  
  #Community-Level Hyperparameters#
  omega ~ dbeta(0.001, 1) #Community inclusion parameter - approximation of Link's scale prior
  
  psi.mean ~ dbeta(1, 1) #Occurence probability intercept - Approximation of Uniform distribution
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  mu.Area ~ dnorm(0,0.1) 
  sd.Area ~ dunif(0, 5)
  tau.Area <- 1/sd.Area^2
  
  mu.Area2 ~ dnorm(0,0.1) 
  sd.Area2 ~ dunif(0, 5)
  tau.Area2 <- 1/sd.Area2^2
  
  p.mean ~ dbeta(1, 1) #Detection probability intercept - Approximation of Uniform distribution
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  rho ~ dunif(-1, 1) 
  tau.eta <- tau.lp/(1 - rho^2)
  
  for(k in 1:M){
    
    #Species-Level Priors#
    w[k] ~ dbern(omega) #Is species k present in the landscape community?
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi) 
    bArea[k] ~ dnorm(mu.Area, tau.Area)
    bArea2[k] ~ dnorm(mu.Area2, tau.Area2)
    mu.eta[k] <- mu.lp + rho * sd.lp/sd.lpsi * (lpsi[k] - mu.lpsi)
    lp[k] ~ dnorm(mu.eta[k], tau.eta)
    logit(p[k]) <- lp[k]
    
    for(i in 1:nPatches){
      
      #Likelihood#
      logit(psi[i, k]) <- lpsi[k] + bArea[k] * Area[i] + bArea2[k] * Area2[i] #Non-conditional occurrence probability 
      psi2[i, k] <- psi[i, k] * w[k]  #Conditional occurrence probability (conditional on inclusion in community)
      zsub[i, k] ~ dbern(psi2[i, k]) #Does species k occurr in the sampled portion of patch i
      y[i, k] ~ dbin(p[k] * zsub[i, k], nSamples) #Model detection records
      psi3[i,k] <- p_comp(p = psi2[i,k], n = n4FullSample[i]) #Probability of species k occuring anywhere in patch i
      z[i, k] ~ dbern(psi3[i,k]) #Predict whether or not the species occurred in the patch by simulating
                                  #sampling of all possible non-overlapping transects within the patch
    }
  }
  
  #Derived Richness Estimates
  for(i in 1:nPatches){
    PatchR_sub[i] <- sum(zsub[i,1:M])
    PatchR_pred[i] <- sum(z[i,1:M])
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
    
    ###Non-Parametric Estimators###
    NPE_pooled <- data$Ypooled #All sampling events pooled
    sample <- sample(1:inputs$Sampling_Events, 1) #Choose 1 sample event at random 
    NPE_single <- data$Y[,, sample] #Isolate single sampling event
    
    NPE_results <- data.frame()
    
    for(j in 1:n_patches){
      
      print(j)
      
      #If SpadeR cannot estimate richness, catch error and return observed richness value
      result_pooled <- tryCatch({
        ChaoSpecies(NPE_pooled[j,], datatype = 'abundance', conf = 0.95)[[3]][c(5,6,9),]
      }, error = function(cond){
        return(matrix(rep(c(sum(ifelse(NPE_pooled[j,] > 0, 1, 0)),NA,NA,NA), 3), nrow = 3, ncol = 4, byrow = T))
      })
      
      result_single <- tryCatch({
        ChaoSpecies(NPE_single[j,], datatype = 'abundance', conf = 0.95)[[3]][c(5,6,9),]
      }, error = function(cond){
        return(matrix(rep(c(sum(ifelse(NPE_single[j,] > 0, 1, 0)),NA,NA,NA), 3), nrow = 3, ncol = 4, byrow = T))
      })
      
      patch_npe <- as.data.frame(rbind(result_pooled, result_single))
      patch_npe$Patch <- j
      patch_npe$Estimator <- rep(c('iChao1', 'ACE', 'Jack2'), 2)
      patch_npe$Samples <- c(rep('Pooled', 3), rep('Single', 3))
      
      names(patch_npe) <- names(NPE_results)
      NPE_results <- rbind(NPE_results, patch_npe)
      
    }
    
    colnames(NPE_results) <- c('Estimate', 'SE.SD', 'LCI', 'UCI', 'Patch', 'Estimator', 'Samples')
    
    
    ###HMSOM###
    #Data preparation
    nSamples <- inputs$Sampling_Events
    Area <- true_data$Patches$Area
    n4FullSample <- round(1/inputs$patch_prop) #Calculate number of non-overlapping transect groups to tesselate whole patches
    
    nSpObs <- data$Observed_Gamma  
    AugSp <- 2 * nSpObs #Augment w/ 2 * no. observed sps
    M <- nSpObs + AugSp
    w <- c(rep(1, nSpObs), rep(NA, AugSp))
    
    y <- data$Yinc 
    yAug <- cbind(y, matrix(0, nrow = n_patches, ncol = AugSp)) #Augment detection array
    z <- ifelse(yAug > 0, 1, NA)
    
    HMSOM_data <- list(y = yAug, zsub = z, z = z, M = M, w = w, n4FullSample = n4FullSample,
                       nPatches = n_patches, nSamples = nSamples, 
                       Area = as.numeric(scale(log(Area))), Area2 = as.numeric(scale(log(Area)^2)))
    
    #Model Fitting
    model <- nimbleModel(HMSOM, constants = HMSOM_data)
    MCMCconf <- configureMCMC(model,
                              monitors = c('PatchR_pred', 'PatchR_sub', 'z', 'zsub'))
    MCMC <- buildMCMC(MCMCconf)
    compModel <- compileNimble(model)
    compMCMC <- compileNimble(MCMC, project = compModel)
    HMSOM_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                         thin = 20, nchains = 4, samplesAsCodaMCMC = T)
    
    ##Extract Posterior richness estimates##
    HMSOM_df <- as.data.frame(as.matrix(HMSOM_out))
    HMSOM_predR <- HMSOM_df[,grepl('PatchR_pred', names(HMSOM_df))]
    HMSOM_subR <- HMSOM_df[,grepl('PatchR_sub', names(HMSOM_df))]
    
    #Predicted richness (whole patches)
    HMSOM_predR_Mn <- apply(HMSOM_predR, 2, mean)
    HMSOM_predR_sd <- apply(HMSOM_predR, 2, sd)
    HMSOM_predR_lci <- apply(HMSOM_predR, 2, quantile, probs = 0.025)
    HMSOM_predR_uci <- apply(HMSOM_predR, 2, quantile, probs = 0.975)
    
    #Richness in sampled areas
    HMSOM_subR_Mn <- apply(HMSOM_subR, 2, mean)
    HMSOM_subR_sd <- apply(HMSOM_subR, 2, sd)
    HMSOM_subR_lci <- apply(HMSOM_subR, 2, quantile, probs = 0.025)
    HMSOM_subR_uci <- apply(HMSOM_subR, 2, quantile, probs = 0.975)
    
    #Collate HMSOM estimates
    HMSOM_predR_results <- as.data.frame(cbind(HMSOM_predR_Mn, HMSOM_predR_sd, HMSOM_predR_lci, HMSOM_predR_uci))
    HMSOM_subR_results <- as.data.frame(cbind(HMSOM_subR_Mn, HMSOM_subR_sd, HMSOM_subR_lci, HMSOM_subR_uci))
    HMSOM_predR_results$Patch <- 
      HMSOM_subR_results$Patch <- 1:n_patches
    HMSOM_predR_results$Estimator <- c(rep('HMSOM_Predicted', n_patches))
    HMSOM_subR_results$Estimator <- c(rep('HMSOM', n_patches))
    HMSOM_predR_results$Samples <- 
      HMSOM_subR_results$Samples <- NA
    colnames(HMSOM_subR_results)[1:4] <- 
      colnames(HMSOM_predR_results)[1:4] <- c('Estimate', 'SE.SD', 'LCI', 'UCI')
    
    ###Observed###
    observed_results <- data.frame(Estimate = rowSums(data$Z > 0),
                                   SE.SD = NA,
                                   LCI = NA,
                                   UCI = NA,
                                   Patch = 1:n_patches, 
                                   Estimator = 'Observed',
                                   Samples = 'Pooled')
    
    #Combine Results
    patchR <- rbind(NPE_results, HMSOM_predR_results, HMSOM_subR_results, observed_results)
    patchR <- merge(true_data$Patches, patchR)
    patchR <- patchR[,c(1:3, 8:9, 4:7)]
    patchR$Estimate_ID <- 1:nrow(patchR)
    patchR$Repeat <- data$Repeat
    
    
    ##################### Beta Diversity ##########################
    
    #Create Pairwise Beta matrices
    comm_single <- data$Y[,,sample] 
    comm_pooled <- data$Ypooled
    
    pairs <- n_patches * (n_patches - 1)/2 #n pairs
    P1 <- P2 <- c()
    
    for(p1 in 1:(n_patches - 1)){ #Create all pairs
      for(p2 in (p1+1):n_patches){
        P1 <- c(P1, p1)
        P2 <- c(P2, p2)
      }
    }
    
    Observed_B <-   
      HMSOM_predB_Ests <-
      HMSOM_subB_Ests <- 
      NPE_B_Ests_Pooled <- 
      NPE_B_Ests_Single <- data.frame(Site_1 = P1,
                                      Site_2 = P2,
                                      Estimate = NA,
                                      SE.SD = NA,
                                      LCI = NA,
                                      UCI = NA)
    
    #True Beta matrix
    true_comm <- true_data$True_PA
    True_B <- data.frame(Site_1 = P1,
                         Site_2 = P2,
                         True_B = NA)
    
    #Extract HMSOM z matrix
    HMSOM_zPred <- HMSOM_df[,grepl('z', names(HMSOM_df))]
    HMSOM_zSub <- HMSOM_df[,grepl('zsub', names(HMSOM_df))]
    
    indices <- vapply(strsplit(colnames(HMSOM_zSub), "[", fixed = T), 
                      `[`, 2, FUN.VALUE = character(1))
    indices <- vapply(strsplit(indices, "]", fixed = T), 
                      `[`, 1, FUN.VALUE = character(1))
    indices <- strsplit(indices, ",")
    ri <- as.numeric(sapply(indices, '[[', 1))
    ci <- as.numeric(sapply(indices, '[[', 2))
    sub_z_mat <- 
      pred_z_mat <- array(NA, dim = c(max(ri), max(ci), nrow(HMSOM_zSub)))
    
    for(ind in 1:length(indices)){ 
      sub_z_mat[ri[ind],ci[ind],] <- HMSOM_zSub[,ind] 
      pred_z_mat[ri[ind],ci[ind],] <- HMSOM_zPred[,ind]
    }
    
    #Estimate beta
    for(k in 1:pairs){
      
      #Observed Beta
      pair_single <- as.data.frame(t(comm_single[c(P1[k],P2[k]),]))
      pair_pooled <- as.data.frame(t(comm_pooled[c(P1[k],P2[k]),]))
      
      Observed_B[k,3] <- sorensen(ifelse(pair_pooled[,1] > 0, 1, 0), 
                                  ifelse(pair_pooled[,2] > 0, 1, 0))
      
      #True beta
      True_B[k,3] <- sorensen(true_comm[P1[k],], true_comm[P2[k],])
      
      #Chao Beta
      #Chao estimators fail when observed shared species == 0
      #catch error and subsequently insert 0 (totally dissimilar) and NAs for variance  
      NPE_B_Ests_Single[k,3:6] <- tryCatch({
        SimilarityPair(pair_single, datatype = 'abundance')$estimated_richness[1,]
      }, error = function(cond){
        return(c(0,NA,NA,NA))
      })
      
      NPE_B_Ests_Pooled[k,3:6] <- tryCatch({
        SimilarityPair(pair_pooled, datatype = 'abundance')$estimated_richness[1,]
      }, error = function(cond){
        return(c(0,NA,NA,NA))
      })
      
      #HMSOM Beta 
      sor_sub_ests <- sor_pred_ests <- vector('numeric', dim(sub_z_mat)[3])
      
      for(iter in 1:dim(sub_z_mat)[3]){ 
        sor_sub_ests[iter] <- sorensen(sub_z_mat[P1[k],,iter], sub_z_mat[P2[k],,iter]) 
        sor_pred_ests[iter] <- sorensen(pred_z_mat[P1[k],,iter], pred_z_mat[P2[k],,iter]) 
      }
      
      #Predicted occurrence (whole patches)
      HMSOM_predB_Ests$Estimate[k] <- mean(sor_pred_ests)
      HMSOM_predB_Ests$SE.SD[k] <- sd(sor_pred_ests)
      HMSOM_predB_Ests$LCI[k] <- quantile(sor_pred_ests, 0.025)
      HMSOM_predB_Ests$UCI[k] <- quantile(sor_pred_ests, 0.975)
      
      #Occurrence in sampled areas
      HMSOM_subB_Ests$Estimate[k] <- mean(sor_sub_ests)
      HMSOM_subB_Ests$SE.SD[k] <- sd(sor_sub_ests)
      HMSOM_subB_Ests$LCI[k] <- quantile(sor_sub_ests, 0.025)
      HMSOM_subB_Ests$UCI[k] <- quantile(sor_sub_ests, 0.975)
      
      cat(k, '\n')
    }
    
    HMSOM_subB_Ests$Estimator <- 'HMSOM'
    HMSOM_predB_Ests$Estimator <- 'HMSOM_Predicted'
    HMSOM_subB_Ests$Samples <- 
      HMSOM_predB_Ests$Samples <- 'All'
    Observed_B$Estimator <- 'Observed'
    Observed_B$Samples <- 'Pooled'
    NPE_B_Ests_Single$Estimator <- 
      NPE_B_Ests_Pooled$Estimator <- 'iChao1'
    NPE_B_Ests_Single$Samples <- 'Single'
    NPE_B_Ests_Pooled$Samples <- 'Pooled'
    
    patchB <- rbind(Observed_B, HMSOM_subB_Ests, HMSOM_predB_Ests, NPE_B_Ests_Pooled, NPE_B_Ests_Single)
    patchB <- merge(patchB, True_B)
    patchB$Estimate_ID <- 1:nrow(patchB)
    patchB$Repeat <- data$Repeat
    
    #Baselga's partitioned beta (SpadeR does not yet provide Simpson dissimilarity so cannot estimate)
    HMSOM_pred_baselga <-
      HMSOM_sub_baselga <- 
      Observed_baselga <- 
      data.frame(Estimator = rep(NA, pairs),
                 Site_1 = P1, Site_2 = P2,
                 Turnover = NA, T_SD = NA, T_UCI = NA, T_LCI = NA,
                 Nestedness = NA, N_SD = NA, N_UCI = NA, N_LCI = NA,
                 True_Turnover = NA, True_Nestedness = NA)
    HMSOM_sub_baselga$Estimator <- 'HMSOM'
    HMSOM_pred_baselga$Estimator <- 'HMSOM_Predicted'
    Observed_baselga$Estimator <- 'Observed'
    
    for(k in 1:pairs){
      
      #True values
      HMSOM_pred_baselga[k, 12:13] <- 
        HMSOM_sub_baselga[k, 12:13] <- 
        Observed_baselga[k, 12:13] <- 
        as.numeric(beta.pair(rbind(true_comm[P1[k],], true_comm[P2[k],]), index.family = 'sorensen')[1:2])
      
      #Observed
      obs_pair <- rbind(comm_pooled[c(P1[k],P2[k]),])
      obs_pair <- ifelse(obs_pair > 0, 1, 0)
      Observed_baselga[k, c(4,8)] <- 
        as.numeric(beta.pair(obs_pair, index.family = 'sorensen')[1:2])
      
      #HMSOM
      pred_baselga_ests <-
        sub_baselga_ests <- matrix(NA, nrow = dim(sub_z_mat)[3], ncol = 2)
      for(iter in 1:dim(sub_z_mat)[3]){ 
        
        sub_baselga_ests[iter,] <- 
          as.numeric(beta.pair(rbind(sub_z_mat[P1[k],,iter], sub_z_mat[P2[k],,iter]), index.family = 'sorensen')[1:2]) 
        
        pred_baselga_ests[iter,] <- 
          as.numeric(beta.pair(rbind(pred_z_mat[P1[k],,iter], pred_z_mat[P2[k],,iter]), index.family = 'sorensen')[1:2]) 
        
      }
      
      #Predicted Baselga beta (whole patch)
      HMSOM_pred_baselga$Turnover[k] <- mean(pred_baselga_ests[,1])
      HMSOM_pred_baselga$T_SD[k] <- sd(pred_baselga_ests[,1])
      HMSOM_pred_baselga$T_UCI[k] <- quantile(pred_baselga_ests[,1], 0.975)
      HMSOM_pred_baselga$T_LCI[k] <- quantile(pred_baselga_ests[,1], 0.025)
      HMSOM_pred_baselga$Nestedness[k] <- mean(pred_baselga_ests[,2])
      HMSOM_pred_baselga$N_SD[k] <- sd(pred_baselga_ests[,2])
      HMSOM_pred_baselga$N_UCI[k] <- quantile(pred_baselga_ests[,2], 0.975)
      HMSOM_pred_baselga$N_LCI[k] <- quantile(pred_baselga_ests[,2], 0.025)
      
      #Baselga beta in sampled areas
      HMSOM_sub_baselga$Turnover[k] <- mean(sub_baselga_ests[,1])
      HMSOM_sub_baselga$T_SD[k] <- sd(sub_baselga_ests[,1])
      HMSOM_sub_baselga$T_UCI[k] <- quantile(sub_baselga_ests[,1], 0.975)
      HMSOM_sub_baselga$T_LCI[k] <- quantile(sub_baselga_ests[,1], 0.025)
      HMSOM_sub_baselga$Nestedness[k] <- mean(sub_baselga_ests[,2])
      HMSOM_sub_baselga$N_SD[k] <- sd(sub_baselga_ests[,2])
      HMSOM_sub_baselga$N_UCI[k] <- quantile(sub_baselga_ests[,2], 0.975)
      HMSOM_sub_baselga$N_LCI[k] <- quantile(sub_baselga_ests[,2], 0.025)
      
      print(k)
      
    }
    
    Baselga_Ests <- rbind(HMSOM_sub_baselga, HMSOM_pred_baselga, Observed_baselga)
    Baselga_Ests$Estimate_ID <- 1:nrow(Baselga_Ests)
    Baselga_Ests$Repeat <- data$Repeat
    
    ######################## SARs ############################
    
    #True SAR
    True_SAR_Data <- unique(patchR[,2:3])
    True_SAR <- lm(log(True_R + 1) ~ log(Area), data = True_SAR_Data)
    True_SAR_Summary <- c(True_SAR$coefficients[1], 
                          True_SAR$coefficients[2])
    
    #Observed SAR
    Observed_SAR_Data <- patchR[patchR$Estimator == 'Observed', c(2,6)]
    Observed_SAR <- lm(log(Estimate + 1) ~ log(Area), data = Observed_SAR_Data)
    Observed_SAR_Summary <- c('Observed', 'LM', Observed_SAR$coefficients[1], 
                              confint(Observed_SAR, '(Intercept)', level = 0.95)[1],
                              confint(Observed_SAR, '(Intercept)', level = 0.95)[2],
                              summary(Observed_SAR)[['coefficients']][1,2],
                              Observed_SAR$coefficients[2],
                              confint(Observed_SAR, 'log(Area)', level = 0.95)[1],
                              confint(Observed_SAR, 'log(Area)', level = 0.95)[2],
                              summary(Observed_SAR)[['coefficients']][2,2])
    
    #Incorporate uncertainty in richness estimates by weighting mean richness by SE/SD
    
    HMSOM_sub_SAR_data <- patchR[patchR$Estimator == 'HMSOM', c(2,6,7)]
    HMSOM_pred_SAR_data <- patchR[patchR$Estimator == 'HMSOM_Predicted', c(2,6,7)]
    
    iChao1_Single_Data <- patchR[patchR$Estimator == 'iChao1' & 
                                   patchR$Samples == 'Single', c(2,6,7)]
    iChao1_Pooled_Data <- patchR[patchR$Estimator == 'iChao1' & 
                                   patchR$Samples == 'Pooled', c(2,6,7)]
    Jack2_Single_Data <- patchR[patchR$Estimator == 'Jack2' & 
                                  patchR$Samples == 'Single', c(2,6,7)]
    Jack2_Pooled_Data <- patchR[patchR$Estimator == 'Jack2' & 
                                  patchR$Samples == 'Pooled', c(2,6,7)]
    ACE_Single_Data <- patchR[patchR$Estimator == 'ACE' & 
                                patchR$Samples == 'Single', c(2,6,7)]
    ACE_Pooled_Data <- patchR[patchR$Estimator == 'ACE' & 
                                patchR$Samples == 'Pooled', c(2,6,7)]
    
    SARs <- list(
      
      #Standard SARs
      HMSOM = lm(log(HMSOM_sub_SAR_data$Estimate + 1) ~ log(HMSOM_sub_SAR_data$Area)),
      HMSOM_Predicted = lm(log(HMSOM_pred_SAR_data$Estimate + 1) ~ log(HMSOM_pred_SAR_data$Area)), 
      iChao1.Pooled = lm(log(iChao1_Single_Data$Estimate + 1) ~ log(iChao1_Single_Data$Area)),
      iChao1.Single = lm(log(iChao1_Pooled_Data$Estimate + 1) ~ log(iChao1_Pooled_Data$Area)),
      Jack2.Pooled = lm(log(Jack2_Single_Data$Estimate + 1) ~ log(Jack2_Single_Data$Area)),
      Jack2.Single = lm(log(Jack2_Pooled_Data$Estimate + 1) ~ log(Jack2_Pooled_Data$Area)),
      ACE.Pooled = lm(log(ACE_Single_Data$Estimate + 1) ~ log(ACE_Single_Data$Area)),
      ACE.Single = lm(log(ACE_Pooled_Data$Estimate + 1) ~ log(ACE_Pooled_Data$Area)),
      
      #Weighted SARs
      HMSOMx = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = HMSOM_sub_SAR_data), #Weight by std dev instead as SE is a frequentist concept
      HMSOM_Predictedx = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = HMSOM_pred_SAR_data),
      iChao1.Pooledx = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = iChao1_Single_Data),
      iChao1.Singlex = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = iChao1_Pooled_Data),
      Jack2.Pooledx = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = Jack2_Single_Data),
      Jack2.Singlex = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = Jack2_Pooled_Data),
      ACE.Pooledx = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = ACE_Single_Data),
      ACE.Singlex = lm(log(Estimate + 1) ~ log(Area), weights = 1/SE.SD, data = ACE_Pooled_Data)
      
    )
    
    
    #Collate results
    SAR_results <- as.data.frame(matrix(nrow = 16, ncol = 10))
    colnames(SAR_results) <- c('Estimator', 'Type', 'c',  'c_LCI', 'c_UCI', 'c_SE',
                               'z', 'z_LCI', 'z_UCI', 'z_SE')
    
    for(n in 1:length(SARs)){
      
      SAR_results$Estimator[n] <- sub('x','', names(SARs)[n])
      SAR_results$c[n] <- as.numeric(coef(SARs[[n]])[1])
      SAR_results$z[n] <- as.numeric(coef(SARs[[n]])[2])
      SAR_results$c_LCI[n] <- confint(SARs[[n]], '(Intercept)', 0.95)[1]
      SAR_results$c_UCI[n] <- confint(SARs[[n]], '(Intercept)', 0.95)[2]
      SAR_results$c_SE[n] <- summary(SARs[[n]])[['coefficients']][1,2]
      SAR_results$z_LCI[n] <- confint(SARs[[n]], 'log(Area)', 0.95)[1]
      SAR_results$z_UCI[n] <- confint(SARs[[n]], 'log(Area)', 0.95)[2]
      SAR_results$z_SE[n] <- summary(SARs[[n]])[['coefficients']][2,2]
      
    }
    
    SAR_results$Type <- c(rep('LM', 8), rep('MA', 8))
    SAR_results <- rbind(Observed_SAR_Summary, SAR_results)
    SAR_results$True_z <-  True_SAR_Summary[2]
    SAR_results$True_c <- True_SAR_Summary[1]
    SAR_results$Estimate_ID <- 1:nrow(SAR_results)
    SAR_results$Repeat <- data$Repeat
    SAR_results[,3:14] <- lapply(SAR_results[,3:14], as.numeric)
    
    #Calculate proportion of patches where NPEs could quantify uncertainty (i.e., observed richness wasn't used instead)
    SAR_results$Non_NA <- NA
    SAR_results$Non_NA[SAR_results$Estimator == 'iChao1_Single'] <- mean(!is.na(iChao1_Single_Data$SE.SD))
    SAR_results$Non_NA[SAR_results$Estimator == 'iChao1_Pooled'] <- mean(!is.na(iChao1_Pooled_Data$SE.SD))
    SAR_results$Non_NA[SAR_results$Estimator == 'Jack2_Single'] <- mean(!is.na(Jack2_Single_Data$SE.SD))
    SAR_results$Non_NA[SAR_results$Estimator == 'Jack2_Pooled'] <- mean(!is.na(Jack2_Pooled_Data$SE.SD))
    SAR_results$Non_NA[SAR_results$Estimator == 'ACE_Single'] <- mean(!is.na(ACE_Single_Data$SE.SD))
    SAR_results$Non_NA[SAR_results$Estimator == 'ACE_Pooled'] <- mean(!is.na(ACE_Pooled_Data$SE.SD))
    
    
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
    BModDat$Estimate[BModDat$Estimate == 1 | BModDat$Estimate > 1] <- 1-1e-3 #Adjust 0 and 1 values for modelling (cannot logit transform 0 or 1)
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
    Bchao_single_dat <- BModDat[BModDat$Estimator == 'iChao1' & 
                                  BModDat$Samples == 'Single',
                                c('Area_Diff', 'Estimate', 'SE.SD')]
    Bchao_pool_dat <- BModDat[BModDat$Estimator == 'iChao1' & 
                                BModDat$Samples == 'Pooled',
                              c('Area_Diff', 'Estimate', 'SE.SD')]
    Bhmsom_sub_dat <- BModDat[BModDat$Estimator == 'HMSOM',
                              c('Area_Diff', 'Estimate', 'SE.SD')]
    Bhmsom_pred_dat <- BModDat[BModDat$Estimator == 'HMSOM_Predicted',
                               c('Area_Diff', 'Estimate', 'SE.SD')]
    
    Bmods <- list(
      
      #Standard models
      iChao1.Single = glm(Estimate ~ Area_Diff, family = binomial, data = Bchao_single_dat),
      iChao1.Pooled = glm(Estimate ~ Area_Diff, family = binomial, data = Bchao_pool_dat),
      HMSOM = glm(Estimate ~ Area_Diff, family = binomial, data = Bhmsom_sub_dat),
      HMSOM_Predicted = glm(Estimate ~ Area_Diff, family = binomial, data = Bhmsom_pred_dat),
      
      #Weighted
      iChao1.Singlex = glm(Estimate ~ Area_Diff, family = binomial, weights = 1/SE.SD, data = Bchao_single_dat),
      iChao1.Pooledx = glm(Estimate ~ Area_Diff, family = binomial, weights = 1/SE.SD, data = Bchao_pool_dat),
      HMSOMx = glm(Estimate ~ Area_Diff, family = binomial, weights = 1/SE.SD, data = Bhmsom_sub_dat),
      HMSOM_Predictedx = glm(Estimate ~ Area_Diff, family = binomial, weights = 1/SE.SD, data = Bhmsom_pred_dat)
      
    )
    
    #Collate results
    BMod_Results <- as.data.frame(matrix(nrow = 8, ncol = 10))
    colnames(BMod_Results) <- c('Estimator', 'Type', 'a',  'a_LCI', 'a_UCI', 'a_SE',
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
    
    BMod_Results$Type <- c(rep('LM', 4), rep('MA', 4))
    BMod_Results <- rbind(BMod_Results, Bobs_mod_summ)
    BMod_Results$True_a <- Btrue_mod_summ[1]
    BMod_Results$True_b <- Btrue_mod_summ[2]
    BMod_Results$Estimate_ID <- 1:nrow(BMod_Results)
    BMod_Results$Repeat <- data$Repeat
    BMod_Results[,3:14] <- lapply(BMod_Results[,3:14], as.numeric)
    
    #Calculate proportion of pairs where iChao could quantify uncertainty (i.e., 0 dissimilarity wasn't assigned instead)
    BMod_Results$Non_NA <- NA
    BMod_Results$Non_NA[BMod_Results$Estimator == 'iChao1_Single'] <- mean(!is.na(Bchao_single_dat$SE.SD))
    BMod_Results$Non_NA[BMod_Results$Estimator == 'iChao1_Pooled'] <- mean(!is.na(Bchao_pool_dat$SE.SD))
    
    
    #################Baselga Models#######################
    #Merge estimates and pairwise area differences
    BasModDat <- merge(Baselga_Ests, dist_mat)
    BasModDat$Turnover[BasModDat$Turnover == 1] <- 1-1e-3
    BasModDat$Nestedness[BasModDat$Nestedness == 1] <- 1-1e-3
    BasModDat$Turnover[BasModDat$Turnover == 0] <- 0+1e-3
    BasModDat$Nestedness[BasModDat$Nestedness == 0] <- 0+1e-3
    BasModDat$True_Turnover[BasModDat$True_Turnover == 1] <- 1-1e-3
    BasModDat$True_Nestedness[BasModDat$True_Nestedness == 1] <- 1-1e-3
    BasModDat$True_Turnover[BasModDat$True_Turnover == 0] <- 0+1e-3
    BasModDat$True_Nestedness[BasModDat$True_Nestedness == 0] <- 0+1e-3
    
    #Extract estimated values
    obs_baselga <- BasModDat[BasModDat$Estimator == 'Observed',]
    hmsom_sub_baselga <- BasModDat[BasModDat$Estimator == 'HMSOM',]
    hmsom_pred_baselga <- BasModDat[BasModDat$Estimator == 'HMSOM_Predicted',]
    
    #True models
    true_nest_mod <- glm(True_Nestedness ~ Area_Diff, family = binomial, data = obs_baselga)
    true_turn_mod <- glm(True_Turnover ~ Area_Diff, family = binomial, data = obs_baselga)
    
    #Observed and HMSOM models
    baselga_mods <- list('Observed-Nestedness.LM' = glm(Nestedness ~ Area_Diff, family = binomial, data = obs_baselga),
                         'Observed-Turnover.LM' = glm(Turnover ~ Area_Diff, family = binomial, data = obs_baselga),
                         'HMSOM-Nestedness.LM' = glm(Nestedness ~ Area_Diff, family = binomial, data = hmsom_sub_baselga),
                         'HMSOM-Turnover.LM' = glm(Turnover ~ Area_Diff, family = binomial, data = hmsom_sub_baselga),
                         'HMSOM_Predicted-Nestedness.LM' = glm(Nestedness ~ Area_Diff, family = binomial, data = hmsom_pred_baselga),
                         'HMSOM_Predicted-Turnover.LM' = glm(Turnover ~ Area_Diff, family = binomial, data = hmsom_pred_baselga),
                         'HMSOM-Nestedness.MA' = glm(Nestedness ~ Area_Diff, family = binomial, weights = 1/N_SD, data = hmsom_sub_baselga),
                         'HMSOM-Turnover.MA' = glm(Turnover ~ Area_Diff, family = binomial, weights = 1/T_SD, data = hmsom_sub_baselga),
                         'HMSOM_Predicted-Nestedness.MA' = glm(Nestedness ~ Area_Diff, family = binomial, weights = 1/N_SD, data = hmsom_pred_baselga),
                         'HMSOM_Predicted-Turnover.MA' = glm(Turnover ~ Area_Diff, family = binomial, weights = 1/T_SD, data = hmsom_pred_baselga))
    
    BaselgaMods_Results <- data.frame(Estimator = rep(NA, 10),
                                      Type = NA, Variable = NA,
                                      a = NA, a_LCI = NA, a_UCI = NA,
                                      b = NA, b_LCI = NA, b_UCI = NA, 
                                      True_a = NA, True_b = NA)
    
    for(n in 1:length(baselga_mods)){
      
      BaselgaMods_Results$Estimator[n] <- sub('-.*', '', names(baselga_mods)[n])
      BaselgaMods_Results$Type[n] <- sub('.*?\\.', '', names(baselga_mods)[n])
      BaselgaMods_Results$Variable[n] <- gsub('.*-(.+)\\..*', '\\1', names(baselga_mods)[n])
      
      BaselgaMods_Results$a[n] <- as.numeric(coef(baselga_mods[[n]])[1])
      BaselgaMods_Results$a_LCI[n] <- confint(baselga_mods[[n]], '(Intercept)', 0.95)[1]
      BaselgaMods_Results$a_UCI[n] <- confint(baselga_mods[[n]], '(Intercept)', 0.95)[2]
      BaselgaMods_Results$b[n] <- as.numeric(coef(baselga_mods[[n]])[2])
      BaselgaMods_Results$b_LCI[n] <- confint(baselga_mods[[n]], 'Area_Diff', 0.95)[1]
      BaselgaMods_Results$b_UCI[n] <-  confint(baselga_mods[[n]], 'Area_Diff', 0.95)[2]
      
    }
    
    BaselgaMods_Results$True_a[BaselgaMods_Results$Variable == 'Nestedness'] <- true_nest_mod$coefficients[1]
    BaselgaMods_Results$True_a[BaselgaMods_Results$Variable == 'Turnover'] <- true_turn_mod$coefficients[1]
    BaselgaMods_Results$True_b[BaselgaMods_Results$Variable == 'Nestedness'] <- true_nest_mod$coefficients[2]
    BaselgaMods_Results$True_b[BaselgaMods_Results$Variable == 'Turnover'] <- true_turn_mod$coefficients[2]
    
    BaselgaMods_Results$Estimate_ID <- 1:nrow(BaselgaMods_Results)
    BaselgaMods_Results$Repeat <- data$Repeat
    
    #####Collate all results and save#####
    results[[i]] <- list(id = id, Parameters = inputs, Repeat = data$Repeat, 
                         Richness_Estimates = patchR, Beta_Estimates = patchB,
                         SAR_Estimates = SAR_results, BetaMod_Estimates = BMod_Results,
                         Baselga_Estimates = Baselga_Ests, BaselgaMod_Estimates = BaselgaMods_Results)
    
    
  }
  
  return(results)
  
}

#####Apply estimators to simulated data#####
#This will be extremely slow but is effectively the same as how the analysis was ran on the UEA HPC cluster, just had 100s more cores on cluster

#Parallel processing
cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)
ntasks <- length(landscape_simulations)
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(x) setTxtProgressBar(pb, x)
opts <- list(progress = progress)

sim_ests <- foreach(i = 1:length(landscape_sims_parsed),
                     .options.snow = opts,
                     .inorder = F,
                     .packages = c('nimble', 'SpadeR', 'metafor', 'betapart', 'reshape2')
  ) %dopar% {
  
  estimate_diversity(landscape_simulations[i])
  
}

stopCluster()
save(sim_ests, file = 'Simulated_Richness_Estimates.RData')
gc()

