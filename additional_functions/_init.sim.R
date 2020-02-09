source(paste(home,"/additional_functions/", "initialization.simulation.R",sep=""))

k_curr <- theta  # Number of sources
initialize_list <- initialization.simulation(params)  

mix_num <- k_curr+1  # Number of mixture components (sources + background)
mu_curr <- initialize_list[[1]]  # Locations of sources (2 x k_curr matrix) - locations initially on edge of image
w <- initialize_list[[2]]  # Relative intensities (vector of length mix_num) 
eparas <- initialize_list[[3]]  # Shape and mean spectral parameters - disregard if spectral_model=="none"
ewt <- initialize_list[[4]]  # Spectral model weights (extended full model) - disregard unless spectral_model=="extended_full"
allocate_curr <- initialize_list[[5]]  # Allocation of photons to sources (and background) - initially obs_num x mix_num matrix of zeros
