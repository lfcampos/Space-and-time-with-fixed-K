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
#  plot(spatial, col = rgb(0, 0, 0, alpha = 2*dpsf), pch = 19)
# hist(energy, prob = TRUE, breaks = 100)
# tmp = seq(min(energy), max(energy), 1)
# lines(tmp, dgamma(tmp, eparasi1[1],eparasi2[1]))
# fix_runs = mcmc_runs
# rjmcmc_run = tt
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
mcmc.spatial.time.energy.full <- function(fix_runs,online_ordering,rjmcmc_run,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin){

  # Number of time breaks is given for all sources
  if(length(num_time_breaks) == 1){
    num_time_breaks = rep(num_time_breaks, k_curr)
  }
  # # checks
  # length(num_time_breaks) == k_curr
  # length(bk) == num_time_breaks + 1
  
  # Standard MCMC updates
  for (t2 in 1:fix_runs){
    # Update photon allocations 
    probs <- matrix(NA,mix_num,obs_num)

    for(i in 1:k_curr){
      dlambda = lambda[[i]][time_bin[[i]]]
      # extract parameters for each photon according to their time bin
      eparas_photon = eparas_all[[i]][time_bin[[i]]]
      eparasi1 = unlist(sapply(eparas_photon, function(x){x[1,1]}))
      eparasi2 = unlist(sapply(eparas_photon, function(x){x[2,1]}))
      # calculate likelihoods
      dpsf = psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,i])
      dE = dgamma(energy,eparasi2,eparasi2/eparasi1)
      probs[i,] = w[i]*dpsf*dE*dlambda
    }

    # update background probabilities
    # subset background lambda accourding to photon time arrival then if enery is bounded above by the max allowed energy.
    dlambda_back = lambda[[mix_num]][time_bin[[mix_num]]][energy <= max_back_energy]

    # prob(background photon) = (mixture weight) * (psf equivalent unif) * (uniform energy dist'n) * (relative time intensity)
    probs[mix_num,energy <= max_back_energy] <- w[mix_num]*(1/img_area)*(1/max_back_energy)*dlambda_back
    probs[mix_num,energy > max_back_energy] <- 0


    allocate_curr <- t(matrix(unlist(lapply(1:obs_num, function(i) rmultinom(1, 1, probs[,i]))),ncol=obs_num))  # Don't need to normalize as rmultinom function does it automatically

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
      eparas_all <- eparas_all[indexmu]
      num_time_breaks[1:k_curr] <- num_time_breaks[indexmu]
      bk[1:k_curr] <- bk[indexmu]
      time_bin[1:k_curr] <- time_bin[indexmu]
    }
    
    # Update weights
    alpha <- rep(wprior,mix_num)
    w <- rdirichlet(1,alpha+count_vector)

    # Update lambda weights
    for(i in 1:mix_num){
      lambda0 <- rep(lambdaprior,num_time_breaks[i])
      count_vector_lambda = table(cut(arrival_time[which(allocate_curr[,i] == 1)], bk[[i]]))
      lambda[[i]] = rdirichlet(1,lambda0 + count_vector_lambda)
    }

    # Update spectral parameters (full model) 
    for (i in 1:k_curr){
      for(k in 1:num_time_breaks[i]){
        index <- which(allocate_curr[,i] == 1 & time_bin[[i]] == k)
        cspectral <- energy[index]
        gmcurr <- eparas_all[[i]][[k]][1,1]
        gacurr <- eparas_all[[i]][[k]][2,1]
        gmprop <- rnorm(1,gmcurr,i*specm_sd)
        gaprop <- rnorm(1,gacurr,speca_sd)
        if ((gmprop > emean.min) & (gmprop < emean.max) & (gaprop > 0)){
          logr <- spectral_post(gmprop,gaprop,NA,NA,NA,cspectral) - spectral_post(gmcurr,gacurr,NA,NA,NA,cspectral)
          u <- runif(1,0,1)
          if (log(u) < logr){
            eparas_all[[i]][[k]][,1] <- c(gmprop,gaprop)
          }
        }
      }
    } # end time break

  }
  
  # Output parameters and log-posterior
  value <- list(c(k_curr,c(mu_curr),c(w)),allocate_curr, eparas_all, ewt_all,lambda, time_bin, bk, num_time_breaks)
  return(value)
}

