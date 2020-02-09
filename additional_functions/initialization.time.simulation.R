initialization.time.simulation <- function(params){
  
  mix_num <- params$num+1
  img_area <- 4*xl*yl
  
  allocate_curr <- matrix(0,obs_num,mix_num)

  for(i in 1:mix_num) allocate_curr[dat[,5] == i, i] = 1
  
  mu_curr <- params$positions
  
  w <- params$N_sim/sum(params$N_sim)

  lambda <- params$lambda.dat

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # Energy Parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  eparas_all <- params$energy_paras

  value <- list(mu_curr,w,eparas_all,ewt_all = NULL,allocate_curr,mix_num,img_area,lambda)
  return(value)
}
  
