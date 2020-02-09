initialization.simulation <- function(params){
  
  mix_num <- params$num+1
  img_area <- 4*xl*yl
  
  allocate_curr <- matrix(0,obs_num,mix_num)

  for(i in 1:mix_num) allocate_curr[dat[,5] == i, i] = 1
  
  mu_curr <- params$positions
  
  w <- params$N_sim/sum(params$N_sim)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # Energy Parameters: for BASCS, start from best fitting gamma
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  require(MASS)
  tmp = tapply(dat[,'eng'], dat[,'source'], function(E){
    m = glm(E~1, family = Gamma(link = "identity"))
    c(m$coefficients[1], gamma.shape(m)[[1]])
  })
  eparas <- do.call('rbind', tmp)[1:2,]

  value <- list(mu_curr,w,eparas,ewt = NULL,allocate_curr,mix_num,img_area)
  return(value)
}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# hist(dat[,'eng'][dat[,'source'] == 1], prob = TRUE, breaks = 100)
# xtic = seq(min(dat[,'eng']), max(dat[,'eng']), 1)
# lines(xtic, dgamma(xtic, eparas[1,2], eparas[1,2]/eparas[1,1]), col = 2)
# hist(dat[,'eng'][dat[,'source'] == 2], prob = TRUE, breaks = 100)
# xtic = seq(min(dat[,'eng']), max(dat[,'eng']), 1)
# lines(xtic, dgamma(xtic, eparas[2,2], eparas[2,2]/eparas[2,1]), col = 2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #