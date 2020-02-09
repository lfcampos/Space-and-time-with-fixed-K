full_log_posterior <- function(mu_curr,allocate_curr,w,k_curr,eparas,ewt,alpha){
  mix_num <- k_curr+1
  loglike <- 0
  like_obs <- matrix(NA,obs_num,mix_num)
  for (j in 1:k_curr){
    like_obs[,j] <- w[j]*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,j])*dgamma(energy,eparas[j,2],eparas[j,2]/eparas[j,1])
  }
  like_obs[,mix_num] <- (w[mix_num]/img_area)*(1/max_back_energy)
  loglike <- sum(log(apply(like_obs,1,sum)))
  components <- matrix(NA,5,1)
  components[1] <- loglike # log likelihood
  components[2] <- -k_curr*log(img_area)   # mu_prior
  components[3] <- log(ddirichlet(w,alpha))   # w_prior
  components[4] <- log(dpois(k_curr,theta))   # k_prior 
  components[5] <- sum(log(dgamma(eparas[,2],ashape,arate)))-k_curr*log(emean.range)# eparas_prior
  value <- c(sum(components),loglike)
  return(value)
}