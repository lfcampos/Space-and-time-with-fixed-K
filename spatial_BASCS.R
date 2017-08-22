# BASCS algorithm: RJMCMC for separating overlapping sources
# SPATIAL DATA ONLY
# David Jones, Harvard University

# See ``Disentangling Overlapping Astronomical Sources using Spatial and Spectral Information"
# Jones, Kashyap, van Dyk (ApJ 2015) for model and method details. 

# Directories
home = "F:/Dropbox/Now/Overlapping Sources/Luis project - time/Disambiguating sources using spatial-time-spectral information/Code/Disambiguating sources 1.5 - space and time with fixed K/" 
results <- paste(home,"results/",sep="")

# Set spectral model type for initialization / function load settings
spectral_model <- "none"

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

# IMPORTANT -- the spatial coordinates have been divided by 200
# so that they are on a scale of 10 arcsecs per unit, which matches our PSF. 
# Using data and a PSF that are on different scales will not give 
# sensible results. 

# Data - spatial and energy
# `spatial': n x 2 matrix of spatial co-ordinates of detected photons
# `energy': n x 1 vector of spectral data (PI channels)
setwd(paste(home,"/data/",sep=""))
load("small_chandra_orion_subset1522.RData")
spatial <- as.matrix(spatial)
colnames(spatial) <- NULL
rownames(spatial) <- NULL

# Basic data plots
par(mfrow=c(1,2))
plot(spatial,xlab="10 arcsecs",ylab="10 arcsecs")
hist(energy,xlab="PI Channel")
dev.off()

# Basic data information
obs_num <- length(spatial[,1])
xlow <- min(spatial[,1])
xup <- max(spatial[,1])
ylow <- min(spatial[,2])
yup <- max(spatial[,2])
yl <- (yup - ylow)/2
xl <- (xup - xlow)/2
img_area <- 4*xl*yl  # Image area - rectangle assumed


##############################################
# Load parametric PSF (I also looked at using images of the PSF but the parametric version is OK for now)
##############################################

# The psf.R function has arguments `pts' (nx2 matrix of spatial points) and 
# `center' (the position at which the PSF is centered). 'psf.norm' the PSF
# normalizing constant that appears in psf.R is calculated automatically 
# below. 

setwd(paste(home,"/psf/",sep=""))
source("psf.R") # 2D King prfile density (King 1962)

function_load_name <- ".R"
slope <- 1.5  # power-law slope
r0 <- 0.6  # core radius (in units of 10 arcsecs)
ellip <- 0.00573533  # ellipticity
off.angle <- 0  # off-axis angle

# Calculate normalizing constant
source("normalize_psf.R")
psf.norm <- 1
psf.norm <- 1/normalize_psf(100)$integral 

# Load C++ version PSF (requires psf.norm calculated above)
sourceCpp('psf_cpp.cpp')


# Visual check PSF that is on the right scale
check <- 0
if (check == 1){
  
  data_density <- kde2d(spatial[,1],spatial[,2],n=round(5*max(diff(range(spatial[,1]))/r0,diff(range(spatial[,1]))/r0)))
  max_index <- which(data_density$z == max(data_density$z),arr.ind=TRUE)
  xmax <- data_density$x[max_index[1]]
  ymax <- data_density$y[max_index[2]]
  
  psf_pos <- c(xmax,ymax)
  i <- 0
  xvec <- seq(xmax-10*r0,xmax+10*r0,0.1)
  yvec <- seq(ymax-10*r0,ymax+10*r0,0.1)
  zvec <- matrix(NA,length(xvec),length(yvec))
  for (xnow in xvec){
    i <- i + 1
    j <- 0
    for (ynow in yvec){
      j <- j + 1
      zvec[i,j] <- psf(c(xnow,ynow),psf_pos)
    }
  }
  plot(spatial,xlab="X",ylab="Y",main="Check scale of PSF against bright source",xlim=range(xvec),ylim=range(yvec))
  contour(xvec,yvec,zvec,col=2,add=TRUE)
  
}

  
####################################################
# Settings
####################################################

# (1) Prior parameter settings
theta <- 14  # NOTE: in MCMC this should be set to the true number of sources
# Prior: K ~ Pois(theta), where K is the numer of sources
# Used in calculation of reported log posterior values (don't change as we'll allow the number of 
# sources to vary eventually)
wprior <- 1  # Prior: w ~ Dirichlet(wprior), where w is a vector giving the relative intensities of the 
# sources (and background)

# (2) MCMC seetings 
adapt_end <- 10000  # Iteration number on which to stop adaptive MCMC location updates
mu_adapt_prop_sd <- 1  # Adaptive MCMC location updates (first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of source i is mu_adapt_prop_var/sqrt(ni)
# where ni is the number of photons currently assigned to 
# source i.                   
mu_fixed_prop_sd <- 0.1  # Location fixed updates (after the first adapt_end iterations):
# the standard deviation of the proposal distribution for 
# updating the location of a source is mu_fixed_prop_var

# (3) Algorithm settings
# Algorithm parameters
burnin_samples <- 1000  # Number of MCMC samples for burnin
mcmc_samples <- 1000  # Number of MCMC samples to obtain
thin <- 1  # How many interations between each stored sample (applies to burnin and main run)
print_interval <- 100  # Number of iterations between print of current iteration number 
# and current number of sources


#############################
# Default initialization
# Users who know the approximate locations of the sources and / or 
# their relative intensities and spectral parameters can input their
# own initialization
#############################

setwd(paste(home,"/additional_functions/",sep=""))
source("initialization.R")

k_curr <- theta  # Number of sources
initialize_list <- initialization(k_curr)  

mix_num <- k_curr+1  # Number of mixture components (sources + background)
mu_curr <- initialize_list[[1]]  # Locations of sources (2 x k_curr matrix) - locations initially on edge of image
w <- initialize_list[[2]]  # Relative intensities (vector of length mix_num) 
allocate_curr <- initialize_list[[5]]  # Allocation of photons to sources (and background) - initially obs_num x mix_num matrix of zeros

plot(spatial,main="Initialization")
points(t(mu_curr),col=2,pch=16)

#############################
# Code required
#############################

setwd(paste(home,"/additional_functions/",sep=""))
# MCMC 
source(paste("mcmc.spatial",function_load_name,sep=""))
source(paste("log_posterior_spatial",function_load_name,sep=""))


#############################
# Initial RJMCMC run
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

initial_run <- mcmc.spatial(burnin_samples,thin,online_ordering,adapt_end,w,allocate_curr,mu_curr,k_curr,mix_num)

# Extract final parameters and allocation matrix
new_paras <- initial_run[[burnin_samples]][[1]]
k_curr <- new_paras[1]
mix_num <- k_curr+1
mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
allocate_curr <- initial_run[[burnin_samples]][[2]]

save(list=ls(),file=paste(results,"burnin_results.RData",sep=""))

plot(spatial,main="After burnin")
points(t(mu_curr),col=2,pch=16)


#############################
# Main RJMCMC run
#############################

online_ordering <- "reference"
no_guess <- k_curr
mu_guess <- matrix(mu_curr[,order(w[-mix_num],decreasing=TRUE)],2,k_curr)
print(mu_guess)
adapt_end <- 0  # Iteration number on which to stop adaptive MCMC location updates

main_run <- mcmc.spatial(mcmc_samples,thin,online_ordering,adapt_end,w,allocate_curr,mu_curr,k_curr,mix_num)

# Extract final parameters and allocation matrix
new_paras <- main_run[[mcmc_samples]][[1]]
k_curr <- new_paras[1]
mix_num <- k_curr+1
mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
allocate_curr <- main_run[[mcmc_samples]][[2]]

# Save output
save(list=ls(),file=paste(results,"main_results.RData",sep=""))

# Show time taken
end <- proc.time()
print(end-start)

# Simple plot of final positions
plot(spatial,main="Final")
points(t(mu_curr),pch=16,col=2)

