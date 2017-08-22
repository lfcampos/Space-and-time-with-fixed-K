# Extended BASCS algorithm: MCMC for separating overlapping sources
# David Jones, Harvard University
# Luis Campos, Harvard University

# See ``Disentangling Overlapping Astronomical Sources using Spatial and Spectral Information"
# Jones, Kashyap, van Dyk (ApJ 2015) for model and method details. 




# args=(commandArgs(TRUE))

# eval(parse(text = args))

# file.name = gsub('./sim_data/', '', file.name)


file.name = 'dat_HBC515_zoom_med.Rda'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# test:
# file.name = 'dat_HBC515_250.Rda'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #



####################################################################################
# FIXED K
# NO RJMCMC VERSION - ignore RJMCMC related notes
####################################################################################

# The allgorithm uses spatial, temporal and spectral information to find the joint posterior distribution
# of the number of sources and their spatial, temporal and spectral parameters.

# RUN TIMES
# The algorithm is computationally intensive. 
# Example: 

############################
# Required R packages
############################

library("MCMCpack")
library("cubature")
library("Rcpp")
library("MASS")

#############################
# Load data
#############################

# Data - spatial and energy
# `spatial': n x 2 matrix of spatial co-ordinates of detected photons
# `energy': n x 1 vector of spectral data (PI channels)
# setwd(paste(home,"/data/",sep=""))
# load("small_chandra_orion_subset1522.RData")

home = getwd()
results <- paste(home,"/results/",sep="")

data_dir = '/data/'
results_dir = '/results/'
initial_run = NULL

out.file.name = gsub('dat', 'analysis', file.name)
load(paste(home, data_dir, file.name, sep = ''))


spatial <- dat[,1:2]
arrival_time <- dat[,4]
energy <- dat[,3]
max_back_energy <- max(energy)


spatial <- as.matrix(spatial)
colnames(spatial) <- NULL
rownames(spatial) <- NULL

# Basic data information
obs_num <- length(spatial[,1])
xlow <- min(spatial[,1])
xup <- max(spatial[,1])
ylow <- min(spatial[,2])
yup <- max(spatial[,2])
yl <- (yup - ylow)/2
xl <- (xup - xlow)/2
img_area <- 4*xl*yl  # Image area - rectangle assumed


###########################
# Options and settings:
###########################

# (1) Time Model: dirichlet/?
#     Spacing: equal/BayesianBlocks
#     num_time_breaks: number of time breaks
time_model <- "dirichlet"
time_spacing <- "equal"
num_time_breaks <- 14

# (1) Spectral model: extended_full / full / none
spectral_model <- "extended_full"

# (2) Load PSF and set parameters 
# The psf.R function has arguments `pts' (nx2 matrix of spatial points) and 
# `center' (the position at which the PSF is centered). 'psf.norm' the PSF
# normalizing constant that appears in psf.R is calculated automatically 
# below. 
setwd(paste(home,"/psf/",sep="")) 

source("psf.R") # 2D King prfile density (King 1962)
# psf.form <- "parametric"
psf.form <- "empirical"

if (psf.form == "parametric"){
  
  function_load_name <- ".R"
  slope <- 1.5  # power-law slope
  r0 <- 0.6  # core radius
  ellip <- 0.00573533  # ellipticity
  off.angle <- 0  # off-axis angle
  
  # Calculate normalizing constant
  source("normalize_psf.R")
  psf.norm <- 1
  psf.norm <- 1/normalize_psf(100)$integral 
  
  # Load C++ version PSF (requires psf.norm calculated above)
  sourceCpp('psf_cpp.cpp')
  
} else {

  # Alternatively load an image of the PSF 
  # User must supply image
  # List of length three:
  # psf_image[[1]] x vector
  # psf_image[[2]] y vector
  # psf_image[[3]] grid of psf values 
  function_load_name <- "_imagepsf.R"
  setwd(paste(home,"/psf/psf_image/",sep=""))
  load("psf_image.RData")
  source("psf_image_function.R") 
  source("adjust_grids.R")
  
  xpsf_grid <- psf_image[[1]] 
  ypsf_grid <- psf_image[[2]] 
  psf_image_values <- psf_image[[3]] 
  
  midpoint_alignment <- adjust_grids()
  xpsf_grid <- midpoint_alignment[[1]]
  ypsf_grid <- midpoint_alignment[[2]]
  psf_image_values <- midpoint_alignment[[3]]  
  sourceCpp('psf_image_function_cpp.cpp')

}

setwd(home)

# (3) Prior parameter settings

lambdaprior <- 1 #Prior:  lambda[[i]] ~ Dirichlet(lambdaprior), where lambda[[i]] is a vector giving the relative 
theta <- 2  # Prior: K ~ Pois(theta), where K is the numer of sources
wprior <- 1  # Prior: w ~ Dirichlet(wprior), where w is a vector giving the relative intensities of the 
# sources (and background)
ashape <- 2  # Prior: a ~ Gamma(ashape,arate), where a is the shape parameter of the spectral model
arate <- 0.5  # See above comment 
ewt_ab <- 2 # ewt ~ Beta(ewt_ab,ewt_ab), where ewt is the spectral weight parameter for the extended full model
emean.min <- min(energy)  # Prior: e ~ Uniform(emean.min,emean.max), where e is the mean parameter of the spectral model
emean.max <- max(energy)  # See above comment
emean.range <- emean.max-emean.min  # The reciporcal of the density of Uniform(emean.min,emean.max)

# (4) MCMC seetings 
adapt_end <- 10000  # Iteration number on which to stop adaptive MCMC location updates
mu_adapt_prop_sd <- 1  # Adaptive MCMC location updates (first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of source i is mu_adapt_prop_var/sqrt(ni)
# where ni is the number of photons currently assigned to 
# source i.                   
mu_fixed_prop_sd <- 0.1  # Location fixed updates (after the first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of a source is mu_fixed_prop_var
specm_sd <- 10  # Proposal: e ~ N(e_current,10*i), where e is the mean parameter of the spectral model, and
# i is the source number (in the inferred sources are ordered so that source 1 is the brightest).
# The proposal is rejected if the proposed value of e is outside (emean.min,emean.max), see (3). 
speca_sd <- 1  # Proposal: a ~ N(a_current,speca_sd), where a is the shape parameter of the spectral model
specwt_sd <- 0.1  # Proposal: ewt ~ inv.logit(rnorm(1,logit(ewt_current),specwt_sd)), where ewt is the weight parameter
# in the extended full model (whcih uses a two gamma spectral model)

# (5) RJMCMC settings

# (6) Algorithm settings
# Algorithm parameters
burnin_mcmc_iters <- 2000  # Number of RJMCMC iterations for the burnin period
main_mcmc_iters <- 25000 # Number of RJMCMC iterations for the main run

# Thinning:
mcmc_runs <- 10  # Number of steps of MCMC to run between each recording or parameters
print_interval <- 100  # Number of iterations between print of current iteration number 

#############################
# Default initialization
# Users who know the approximate locations of the sources and / or 
# their relative intensities and spectral parameters can input their
# own initialization
#############################

source(paste(home,"/additional_functions/", "initialization.time.R",sep=""))

k_curr <- theta  # Number of sources
initialize_list <- initialization.time(k_curr, num_time_breaks)  

mix_num <- k_curr+1  # Number of mixture components (sources + background)
mu_curr <- initialize_list[[1]]  # Locations of sources (2 x k_curr matrix) - locations initially on edge of image
w <- initialize_list[[2]]  # Relative intensities (vector of length mix_num) 
eparas_all <- initialize_list[[3]]  # Shape and mean spectral parameters - disregard if spectral_model=="none"
ewt_all <- initialize_list[[4]]  # Spectral model weights (extended full model) - disregard unless spectral_model=="extended_full"
allocate_curr <- initialize_list[[5]]  # Allocation of photons to sources (and background) - initially obs_num x mix_num matrix of zeros
lambda <- initialize_list[[8]]  # Relative Intensities of time arrival distribution


bk = list()
time_bin = list()
if(time_spacing == 'equal'){
  for(i in 1:(k_curr + 1)){
    bk[[i]] = seq(min(arrival_time), max(arrival_time), length.out = num_time_breaks + 1)
    time_bin[[i]] = cut(arrival_time, bk[[i]], include.lowest = TRUE)
    levels(time_bin[[i]]) = 1:num_time_breaks
    time_bin[[i]] = as.numeric(time_bin[[i]])
  }
}




#############################
# Code required
#############################

setwd(paste(home, "/additional_functions/", sep=""))
# MCMC 
source(paste("mcmc.extended.full",function_load_name,sep=""))
source(paste("mcmc.spatial.time.energy.extended.full",function_load_name,sep=""))
# source(paste("mcmc.full",function_load_name,sep=""))
source(paste("extended_full_log_posterior",function_load_name,sep=""))
# source(paste("full_log_posterior",function_load_name,sep=""))
source("logit.R")
source("inv.logit.R")
if (spectral_model=="full"){
  fixedk_mcmc <- mcmc.full
  source("spectral_post_full.R")
} else {
  fixedk_mcmc <- mcmc.extended.full
  source("spectral_post_extended_full.R")
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Luis: fix the likelihood calculation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setwd(home)

#############################
# Initial MCMC run
#############################

online_ordering <- "none"  # How the sources should be ordered after each iteration. Options: "none" / "reference"
# The option "reference" means that the sources will be ordered by their suspected
# relative intensities. In particular, the current sources will be matched 
# (as near as possible) to reference positions. The reference positions are typically
# obtained by running the RJMCMC with no ordering until convergence and then using 
# the final iteration as the reference for new draws. A reference works well if the source
# position posteriors do not have substantial overlap. 

# Record start time
start <- proc.time()

adapt_end <- round(burnin_mcmc_iters/2)
paras <- list()
store_log_post <- list()
store_log_like <- list()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Initial run: run BASCS on entire dataset to get a good initial allocation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
eparas = eparas_all[[1]]
ewt = ewt_all[[1]]

# burnin_mcmc_iters = 5000

for (tt in 1:burnin_mcmc_iters){
  
  # Print iteration at intervals
  if (tt/print_interval == round(tt/print_interval)){
    print(paste("Iteration number: ",tt,sep=""))
    print(paste("Number of sources: ",k_curr,sep=""))
  }
  initial_run <- mcmc.extended.full(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)
  
  # Extract parameters and allocation matrix after burnin
  new_paras <- initial_run[[1]]
  k_curr <- new_paras[1]
  mix_num <- k_curr+1
  mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
  w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
  eparas <- matrix(new_paras[(3*k_curr+3):(7*k_curr+2)],nrow=k_curr)  
  ewt <- new_paras[(7*k_curr+3):(8*k_curr+2)]
  allocate_curr <- initial_run[[2]]
  
  paras[[tt]] <- new_paras
  store_log_post[[tt]] <- initial_run[[3]][1]
  store_log_like[[tt]] <- initial_run[[3]][2]
  
}



initial_run <-  list(paras,allocate_curr,store_log_post,store_log_like)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Initial run: Run Extended BASCS 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

for (tt in 1:burnin_mcmc_iters){
  # Print iteration at intervals
  if (tt/print_interval == round(tt/print_interval)){
    print(paste("Iteration number: ",tt,sep=""))
    # print(paste("Number of sources: ",k_curr,sep=""))
    if(file.name == 'dat_HBC515.Rda') try(save(initial_run_time, file = './results/initial_run_time_temp.Rda'))
  }
  initial_run_time <- mcmc.spatial.time.energy.extended.full(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks)
  
  # Extract parameters and allocation matrix after burnin
  new_paras <- initial_run_time[[1]]
  k_curr <- new_paras[1]
  mix_num <- k_curr+1
  mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
  w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]

  eparas_all <- initial_run_time[[3]] 
  ewt_all <- initial_run_time[[4]] 

  allocate_curr <- initial_run_time[[2]]
  
  paras[[tt]] <- list(par = new_paras, eparas_all = eparas_all, ewt_all = ewt_all)

}

initial_run_time <-  list(paras,allocate_curr)



#############################
# Main MCMC run
#############################

online_ordering <- "reference"
no_guess <- k_curr
mu_guess <- mu_curr[,order(w[-mix_num],decreasing=TRUE)]
adapt_end <- 0  # Iteration number on which to stop adaptive MCMC location updates

paras <- list()
store_log_post <- list()
store_log_like <- list()
store_alloc <- list()

for (tt in 1:main_mcmc_iters){
  
  # Print iteration at intervals
  if (tt/print_interval == round(tt/print_interval)){
    print(paste("Iteration number: ",tt,sep=""))
    # print(paste("Number of sources: ",k_curr,sep=""))
    try(save(main_run_time, file = './results/main_run_time_temp.Rda'))
  }
  
  main_run_time <- mcmc.spatial.time.energy.extended.full(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks)

  # Extract parameters and allocation matrix after burnin
  new_paras <- main_run_time[[1]]
  k_curr <- new_paras[1]
  mix_num <- k_curr+1
  mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
  w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]

  eparas_all <- main_run_time[[3]] 
  ewt_all <- main_run_time[[4]] 

  allocate_curr <- main_run_time[[2]]
  # ouput allocation for complete data only
  if(TRUE){
    save(allocate_curr,file=paste(results, 'alloc/', 'alloc', tt, '.Rda',sep=""))
  } 

  paras[[tt]] <- list(par = new_paras, eparas_all = eparas_all, ewt_all = ewt_all)

}
main_run_time <-  list(paras,allocate_curr)

# Save output
save(list=ls(),file=paste(results,out.file.name,sep=""))

# Show time taken
end <- proc.time()
print(end-start)

# # Simple plot of final positions
# plot(spatial)
# points(mu_curr,pch=16,col=2)



# lapply(paras, function(x) points(rbind(x[c(2, 4)], x[c(3, 5)]), pch = 19, col = 2, cex = 0.3))

# plot(spatial, pch = 19, col = allocate_curr %*% c(0, 1, 0), cex = 0.3)
# points(spatial, pch = 19, col = allocate_curr %*% c(2, 0, 0), cex = 0.3)



# plot(spatial, pch = 19, cex = 0.3, col = allocate_curr %*% c(0, 1, 0))
# points(t(mu_curr),pch=16,col=2)
