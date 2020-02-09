# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# input:
#   w: source intensity (current) of length k_curr + 1
#   allocate_curr: vector of allocations (current) of length obs_num
#   mu_curr: matrix of source locations (current)
#   eparas_all: list of energy parameters of length number of breakpoints per source
#   ewt_all: list of energy parameters of length number of breakpoints per source
#   bk: list of breakpoints (one vec of ln num_time_breaks + 1 for each souce)
#   num_time_breaks: number of breaks considered
#   lambda: current rate vectors for time arrival dirichlet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Testing: plot(energy, (ewti*dgamma(energy,eparasi3,eparasi3/eparasi1)+(1-ewti)*dgamma(energy,eparasi4,eparasi4/eparasi2)), col = time_bin[[i]], pch = 19,  cex = 0.2)
# fix_runs = mcmc_runs
# rjmcmc_run = tt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
mcmc.spatial.time.energy.extended.full <- function(fix_runs,online_ordering,rjmcmc_run,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks, lambda){
  
  # Standard MCMC updates
  for (t2 in 1:fix_runs){
    
    # Update photon allocations 
    probs <- matrix(NA,mix_num,obs_num)

    for(i in 1:k_curr){
      dlambda = lambda[[i]][time_bin[[i]]]
      # extract parameters for each photon according to their time bin
      eparas_photon = eparas_all[time_bin[[i]]]
      ewt_photon = ewt_all[time_bin[[i]]]
      eparasi3 = unlist(sapply(eparas_photon, function(x){x[i,3]}))
      eparasi1 = unlist(sapply(eparas_photon, function(x){x[i,1]}))
      eparasi4 = unlist(sapply(eparas_photon, function(x){x[i,4]}))
      eparasi2 = unlist(sapply(eparas_photon, function(x){x[i,2]}))
      ewti = unlist(sapply(ewt_photon, function(x){x[i]}))
      # calculate likelihoods
      dpsf = psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,i])
      dE = (ewti*dgamma(energy,eparasi3,eparasi3/eparasi1)+(1-ewti)*dgamma(energy,eparasi4,eparasi4/eparasi2))
      probs[i,] = w[i]*dpsf*dE*dlambda
    }

    # update background probabilities
    # subset background lambda accourding to photon time arrival then if enery is bounded above by the max allowed energy.
    dlambda_back = lambda[[mix_num]][time_bin[[mix_num]]][energy <= max_back_energy]

    # prob(background photon) = (mixture weight) * (psf equivalent unif) * (uniform energy dist'n) * (relative time intensity)
    probs[mix_num,energy <= max_back_energy] <- w[mix_num]*(1/img_area)*(1/max_back_energy)*dlambda_back
    probs[mix_num,energy > max_back_energy] <- 0



    allocate_curr <- t(matrix(unlist(lapply(1:obs_num,function(i) rmultinom(1, 1, probs[,i]))),ncol=obs_num))  # Don't need to normalize as rmultinom function does it automatically

    # Counts vector
    mix_num <- k_curr+1
    count_vector <- matrix(NA,mix_num,1)
    count_vector[1:mix_num] <- apply(allocate_curr[,1:mix_num],2,sum)
    
    # Update positions
    mu_prop <- mu_curr 
    for (i in 1:k_curr){      
      index <- allocate_curr[,i]==1
      if (count_vector[i]>0){
        # Adaptive version (eventually ended to ensure convegence)
        if (rjmcmc_run < adapt_end){
          mu_prop[,i] <- rnorm(2,mu_curr[,i],mu_adapt_prop_sd/sqrt(count_vector[i]))
        }
        # Non-adaptive version
        if (rjmcmc_run >= adapt_end){
          mu_prop[,i] <- rnorm(2,mu_curr[,i],mu_fixed_prop_sd)
        }
        logr <- sum(log(psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[index,],mu_prop[,i])))-sum(log(psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial[index,],mu_curr[,i])))
        u <- runif(1,0,1)
        if(is.na(logr)==0){
          if (log(u) < logr){
            mu_curr[,i] <- mu_prop[,i]
          }
        }
      }
      # Have to make sure that sources without phoons assigned move (doesn't effect likelihood) 
      if (count_vector[i]==0){
        if (rjmcmc_run < adapt_end){
          mu_curr[,i] <- c(runif(1,xlow,xup),runif(1,ylow,yup))
        } else {
          mu_curr[,i] <- rnorm(2,mu_curr[,i],mu_fixed_prop_sd)
        }
      }
    }
    
    # Order parameters by suspected sources intensities (associated with particular positions)
    if (k_curr > 1 & online_ordering =="reference"){
      to_order <- min(no_guess,k_curr)
      next_index <- which.min(apply((mu_curr-mu_guess[,1])^2,2,sum))
      next_index_store <- next_index
      for (i in 2:to_order){
        next_order <- order(apply((mu_curr-mu_guess[,i])^2,2,sum))
        next_index_all <- setdiff(next_order,next_index_store)
        next_index_store <- c(next_index_store,next_index_all[1])
      }
      indexmu <- c(next_index_store,setdiff(1:k_curr,next_index_store))
      mu_curr <- mu_curr[,indexmu]
      count_vector[1:k_curr] <- count_vector[indexmu]
      allocate_curr[,1:k_curr] <- allocate_curr[,indexmu]
      for(k in 1:num_time_breaks){
        eparas_all[[k]] <- eparas_all[[k]][indexmu,] 
        ewt_all[[k]] <- ewt_all[[k]][indexmu]
      }
    }
    
    # Update weights
    alpha <- rep(wprior,mix_num)
    w <- rdirichlet(1,alpha+count_vector)

    # Update lambda weights
    for(i in 1:mix_num){
      lambda0 <- rep(lambdaprior,num_time_breaks)
      count_vector_lambda = table(cut(arrival_time[which(allocate_curr[,i] == 1)], bk[[i]]))
      lambda[[i]] = rdirichlet(1,lambda0 + count_vector_lambda)
    }

    # Update spectral parameters (extended full model) 
    for(k in 1:num_time_breaks){
      for (i in 1:k_curr){
        index <- which(allocate_curr[,i] == 1 & time_bin[[i]] == k)
        cspatial <- energy[index]
        gm1curr <- eparas_all[[k]][i,1]
        gm2curr <- eparas_all[[k]][i,2]
        ga1curr <- eparas_all[[k]][i,3]
        ga2curr <- eparas_all[[k]][i,4]
        ewtcurr <- ewt_all[[k]][i]
        ewtprop <- inv.logit(rnorm(1,logit(ewtcurr),specwt_sd))
        gm1prop <- rnorm(1,gm1curr,i*specm_sd)
        ga1prop <- rnorm(1,ga1curr,speca_sd)
        if ((gm1prop > emean.min) & (gm1prop < emean.max) & (ga1prop > 0)){
          logr <- spectral_post(gm1prop,ga1prop,gm2curr,ga2curr,ewtprop,cspatial) - spectral_post(gm1curr,ga1curr,gm2curr,ga2curr,ewtcurr,cspatial)
          u <- runif(1,0,1)
          if (log(u) < logr){
            eparas_all[[k]][i,c(1,3)] <- c(gm1prop,ga1prop)
            ewt_all[[k]][i] <- ewtprop
            gm1curr <- eparas_all[[k]][i,1]
            ga1curr <- eparas_all[[k]][i,3]
            ewtcurr <- ewt_all[[k]][i]
          }
        }
        gm2prop <- rnorm(1,gm2curr,i*specm_sd)
        ga2prop <- rnorm(1,ga2curr,speca_sd)
        if ((gm2prop > emean.min) & (gm2prop < emean.max) & (ga2prop > 0)){
          logr <- spectral_post(gm1curr,ga1curr,gm2prop,ga2prop,ewtcurr,cspatial) - spectral_post(gm1curr,ga1curr,gm2curr,ga2curr,ewtcurr,cspatial)
          u <- runif(1,0,1)
          if (log(u) < logr){
            eparas_all[[k]][i,c(2,4)] <- c(gm2prop,ga2prop)
          }
        }
      }

      # Order Gammas for identifiability and combine conditions
      epara_order <- apply(eparas_all[[k]][,c(1,2)],1,order)
      epara_order_all <- rbind(epara_order,epara_order+2)
      for (i in 1:k_curr){
        eparas_all[[k]][i,] <- eparas_all[[k]][i,epara_order_all[,i]] 
      }
      for (i in 1:k_curr){
        if (sum(epara_order[,i] != c(1,2))>0){
          ewt_all[[k]][i] <- 1- ewt_all[[k]][i]
        }
      }

    } # end time break

  }
  
  # Output parameters and log-posterior
  value <- list(c(k_curr,c(mu_curr),c(w)),allocate_curr, eparas_all, ewt_all,lambda)
  return(value)
}
