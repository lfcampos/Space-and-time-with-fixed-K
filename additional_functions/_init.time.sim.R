source(paste(home,"/additional_functions/", "initialization.time.simulation.R",sep=""))

k_curr <- theta  # Number of sources
initialize_list <- initialization.time.simulation(params)  

mix_num <- k_curr+1  # Number of mixture components (sources + background)
mu_curr <- t(initialize_list[[1]])  # Locations of sources (2 x k_curr matrix) - locations initially on edge of image
w <- initialize_list[[2]]  # Relative intensities (vector of length mix_num) 
eparas_all <- initialize_list[[3]]  # Shape and mean spectral parameters - disregard if spectral_model=="none"
ewt_all <- initialize_list[[4]]  # Spectral model weights (extended full model) - disregard unless spectral_model=="extended_full"
allocate_curr <- initialize_list[[5]]  # Allocation of photons to sources (and background) - initially obs_num x mix_num matrix of zeros
lambda <- initialize_list[[8]]  # Relative Intensities of time arrival distribution


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Initialize breakpoints and time_bins
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

bk = list()
time_bin = list()

if(time_spacing == 'test'){
  q = as.numeric(q)

  num_time_breaks = c(sapply(params$breakpoints, length)[1:2]-1, 1)

  bk =  params$breakpoints
  bk[[3]] = c(min(params$breakpoints[[3]]), max(params$breakpoints[[3]]))

  # add jitter to the middle breakpoints
  for(i in 1:(k_curr)){
    # bk[[i]][-c(1, num_time_breaks[i] + 1)] = bk[[i]][-c(1, num_time_breaks[i] + 1)] + rnorm(num_time_breaks[i]-1, 0, q)
    bk[[i]][-c(1, num_time_breaks[i] + 1)] = bk[[i]][-c(1, num_time_breaks[i] + 1)] + q 
  }
  for(i in 1:(k_curr + 1)){
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks[i]
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }

  eparas_all_tmp = list()
  lambda_tmp = list()

  for(i in 1:(k_curr)){
    eparas_all_tmp[[i]] = list()
    for(k in 1:num_time_breaks[i]){
      eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
    }
  }
  eparas_all = eparas_all_tmp
  # initialize lambda
  for(i in 1:(k_curr+1)){
    lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  }

  eparas_all = eparas_all_tmp
  lambda = lambda_tmp

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Merging Breakpoints: we test merging a few breakpoints
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if(time_spacing == 'test_merging'){
  # test case have 5 bins, let's replicate what BB is doing
  # [1] 0.0006081108 0.2022975916 0.6218380892 0.9987439354
  # [1] 0.0006081108 0.3629292103 0.9987439354
  # [1] 0.0006081108 0.9987439354

  # bright: merge 

  num_time_breaks = c(3, 2, 1)

  bk =  params$breakpoints
  bk[[1]] = bk[[1]][c(1, 2, 4, 6)]
  bk[[2]] = bk[[2]][c(1, 3, 6)]
  bk[[3]] = c(min(params$breakpoints[[3]]), max(params$breakpoints[[3]]))

  for(i in 1:(k_curr + 1)){
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks[i]
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }

  eparas_all_tmp = list()
  lambda_tmp = list()

  for(i in 1:(k_curr)){
    eparas_all_tmp[[i]] = list()
    for(k in 1:num_time_breaks[i]){
      eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
    }
  }
  eparas_all = eparas_all_tmp
  # initialize lambda
  for(i in 1:(k_curr+1)){
    lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  }

  eparas_all = eparas_all_tmp
  lambda = lambda_tmp

}

if(time_spacing == 'test_merging_one'){
  # test case have 5 bins, let's merge just one block see what effect is

  num_time_breaks = c(4, 4, 1)

  bk =  params$breakpoints
  bk[[1]] = bk[[1]][-2]
  bk[[2]] = bk[[2]][-5]
  bk[[3]] = c(min(params$breakpoints[[3]]), max(params$breakpoints[[3]]))

  for(i in 1:(k_curr + 1)){
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks[i]
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }

  eparas_all_tmp = list()
  lambda_tmp = list()

  for(i in 1:(k_curr)){
    eparas_all_tmp[[i]] = list()
    for(k in 1:num_time_breaks[i]){
      eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
    }
  }
  eparas_all = eparas_all_tmp
  # initialize lambda
  for(i in 1:(k_curr+1)){
    lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  }

  eparas_all = eparas_all_tmp
  lambda = lambda_tmp

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if(time_spacing == 'simulation'){

  num_time_breaks = c(sapply(params$breakpoints, length)[1:2]-1, 1)

  bk =  params$breakpoints
  bk[[3]] = c(min(params$breakpoints[[3]]), max(params$breakpoints[[3]]))

  for(i in 1:(k_curr + 1)){
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks[i]
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }

  eparas_all_tmp = list()
  lambda_tmp = list()

  for(i in 1:(k_curr)){
    eparas_all_tmp[[i]] = list()
    for(k in 1:num_time_breaks[i]){
      eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
    }
  }
  eparas_all = eparas_all_tmp
  # initialize lambda
  for(i in 1:(k_curr+1)){
    lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  }

  eparas_all = eparas_all_tmp
  lambda = lambda_tmp

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Equal time spacing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if(time_spacing == 'equal'){

  for(i in 1:(k_curr + 1)){
    bk[[i]] = seq(min(arrival_time), max(arrival_time), length.out = num_time_breaks + 1)
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }
  num_time_breaks = rep(num_time_breaks, k_curr + 1)

  # rearrange initial (alpha, gamma)
  eparas_all_tmp = list()
  lambda_tmp = list()
  for(i in 1:(k_curr)){
    eparas_all_tmp[[i]] = list()
    for(k in 1:num_time_breaks[i]){
      eparas_all_tmp[[i]][[k]] = matrix(eparas_all[[k]][i,], nrow = 2, ncol = 1)
    }
  }
  eparas_all = eparas_all_tmp
  # initialize lambda
  for(i in 1:(k_curr+1)){
    lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  }

  eparas_all = eparas_all_tmp
  lambda = lambda_tmp

}

if(time_spacing == 'bayesian.blocks'){
  library(reticulate)
  astropy <- import("astropy")

  num_time_breaks = rep(NA, k_curr + 1)
  for(i in 1:(k_curr)){
    bk[[i]] = astropy$stats$bayesian_blocks(c(min(arrival_time), arrival_time[allocate_curr[,i]==1], max(arrival_time)))
    num_time_breaks[i] = length(bk[[i]]) - 1

    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:length(levels(time_bin[[i]]))
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }
  # background photons
  bk[[k_curr + 1]] = c(min(arrival_time), max(arrival_time))
  num_time_breaks[k_curr + 1] = 1
  time_bin[[k_curr + 1]] = rep(1, length(arrival_time))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  # Reformat other initialization parameters based on these breakpoints
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  eparas_all_tmp = list()
  lambda_tmp = list()
  # gamma MLE
  f = function(theta){  
    alpha = theta[1]
    gamma = theta[2]
    -sum(dgamma(eik, gamma, gamma/alpha, log = TRUE))
  }

  # initialize (alpha, gamma) at MLE for each source/timepoint
  for(i in 1:(k_curr)){
    eparas_all_tmp[[i]] = list()
    for(k in 1:num_time_breaks[i]){
      eik = energy[allocate_curr[,i]==1 & time_bin[[i]] == k]
      if(length(eik) < 20) eik = energy[allocate_curr[,i]==1]
      MLE = optim(eparas_all[[i]][1,], f)$par
      eparas_all_tmp[[i]][[k]] = matrix(MLE, nrow = 2, ncol = 1)
    }
  }
  # initialize lambda
  for(i in 1:(k_curr+1)){
    lambda_tmp[[i]] = prop.table(table(time_bin[[i]][allocate_curr[,i]==1]))
  }

  eparas_all = eparas_all_tmp
  lambda = lambda_tmp

}


