initialization.time <- function(k_curr, num_time_breaks){
  
  mix_num <- k_curr+1
  img_area <- 4*xl*yl
  
  allocate_curr <- matrix(0,obs_num,mix_num)
  if (k_curr > 1){
    if (k_curr %% 2 == 0){
      first <- k_curr/2
    }
    if (k_curr %% 2 == 1){
      first <- (k_curr+1)/2
    }
    mu_curr <- matrix(cbind(rbind(seq(xlow,xup,length.out=first),ylow),rbind(seq(xlow,xup,length.out=k_curr-first),yup)),nrow=2,ncol=k_curr)
  }
  if (k_curr == 1){
    xposs <- c(xup,xlow)
    yposs <- c(yup,ylow)
    xfirst <- xposs[sample(1:2,1)]
    yfirst <- yposs[sample(1:2,1)]
    mu_curr <- matrix(c(xfirst,yfirst),nrow=2,ncol=k_curr)
  }
  
  w <- rep(1/mix_num,mix_num)

  lambda <- replicate(mix_num, {rep(1/num_time_breaks,num_time_breaks)}, simplify = FALSE)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # Energy Parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  energy_parameters_all = replicate(num_time_breaks, {
    ewt <- rep(0.5,k_curr)
    no.eparas <- ifelse(spectral_model == "full",2,4)
    if (spectral_model != "none"){
      eparas <- matrix(NA,k_curr,no.eparas)
      for (i in 1:k_curr){
        gmeans <- runif(2,emean.min,emean.max)
        if (spectral_model == "full"){
          eparas[i,1] <- gmeans[1]
        } else {
          eparas[i,1] <- min(gmeans)
          eparas[i,2] <- max(gmeans)
        }
      }
      if (spectral_model == "full"){
        eparas[,2] <- 5
      } else {
        eparas[,3:4] <- 5
      }
    }
    list(eparas = eparas, ewt = ewt)
  }, simplify = FALSE)
  allocate_curr <- matrix(0,obs_num,mix_num)
  
  eparas_all = lapply(energy_parameters_all, function(x) x[[1]])
  ewt_all = lapply(energy_parameters_all, function(x) x[[2]])

  value <- list(mu_curr,w,eparas_all,ewt_all,allocate_curr,mix_num,img_area,lambda)
  return(value)
}
  
