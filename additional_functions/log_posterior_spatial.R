log_posterior_spatial <- function(mu_curr,allocate_curr,w,k_curr){
  mix_num <- k_curr+1
  loglike <- 0
  like_obs <- matrix(NA,obs_num,mix_num)
  for (j in 1:k_curr){
    like_obs[,j] <- w[j]*psf_cpp(off.angle,ellip,slope,psf.norm,r0,spatial,mu_curr[,j])
  }
  like_obs[,mix_num] <- (w[mix_num]/img_area)
  loglike <- sum(log(apply(like_obs,1,sum)))
  components <- matrix(NA,4,1)
  components[1] <- loglike # log likelihood
  components[2] <- -k_curr*log(img_area)   # mu_prior
  alpha <- rep(wprior,mix_num)
  components[3] <- log(ddirichlet(w,alpha))   # w_prior
  components[4] <- log(dpois(k_curr,theta))   # k_prior 
  value <- c(sum(components),loglike)
  return(value)
}