###########################
# Options and settings:
###########################

# (0) Time Model: dirichlet
#     Spacing: equal, bayesian.blocks
#     num_time_breaks: number of time breaks (for "equal")
time_model <- "dirichlet"

time_spacing <- "equal"
num_time_breaks <- 5
# time_spacing <- "bayesian.blocks"
# time_spacing <- "simulation"
# time_spacing <- "test"    # we shift the middle blocks by q
# time_spacing <- "test_merging_one"  # we merge just two blocks
# time_spacing <- "test_merging"  # we merge some blocks (similar to what bb does, needs 5 blocks per source)



# (1) Spectral model: extended_full / full / none
# Note, currenly extended full does not handle different number of time bins per source
spectral_model <- "full"


# (2) Set psf parameters 
# The psf.R function has arguments `pts' (nx2 matrix of spatial points) and 
# `center' (the position at which the PSF is centered). 'psf.norm' the PSF
# normalizing constant that appears in psf.R is calculated automatically 
# below. 
psf.form <- "parametric"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setwd(paste(home,"/psf/",sep="")) 
source("psf.R") # 2D King prfile density (King 1962)

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
  
}

setwd(home)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Load required functions MCMC and Posterior Calculations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setwd(paste(home, "/additional_functions/", sep=""))
source(paste("mcmc.full",function_load_name,sep=""))
source(paste("mcmc.spatial.time.energy.full",function_load_name,sep=""))
source(paste("mcmc.extended.full",function_load_name,sep=""))
source(paste("mcmc.spatial.time.energy.extended.full",function_load_name,sep=""))
source(paste("extended_full_log_posterior",function_load_name,sep=""))
source(paste("full_log_posterior",function_load_name,sep=""))
source("logit.R")
source("inv.logit.R")

if (spectral_model=="full"){
  fixedk_mcmc <- mcmc.full
  fixedk_eBASCS_mcmc <- mcmc.spatial.time.energy.full
  source("spectral_post_full.R")
} else {
  fixedk_mcmc <- mcmc.extended.full
  fixedk_eBASCS_mcmc <- mcmc.spatial.time.energy.extended.full
  source("spectral_post_extended_full.R")
}

setwd(home)