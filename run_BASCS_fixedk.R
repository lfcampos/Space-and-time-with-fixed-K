# BASCS algorithm: RJMCMC for separating overlapping sources
# David Jones, Harvard University

# See ``Disentangling Overlapping Astronomical Sources using Spatial and Spectral Information"
# Jones, Kashyap, van Dyk (ApJ 2015) for model and method details. 




args=(commandArgs(TRUE))
eval(parse(text = args))

file.name = gsub('./sim_data/', '', file.name)

start <- proc.time()

####################################################################################
# FIXED K
# NO RJMCMC VERSION - ignore RJMCMC related notes
####################################################################################

############################
# Required R packages
############################

library("MASS")
library("MCMCpack")
library("cubature")
library("Rcpp")

#############################
# Load data
#############################

home = getwd()
results <- paste(home,"/results/",sep="")
data_dir = '/sim_data/'

initial_run = NULL

out.file.name = gsub('simulation', 'analysis', file.name)
load(paste(home, data_dir, file.name, sep = ''))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (1) Split data into component pieces and calculate image parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_clean_data.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (2) Set model options, number of breaks, psf, etc.
# Note: some options set for eBASCS, but ignored in this MCMC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_model_options.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (3) Prior parameter settings
# Note: some options set for eBASCS, but ignored in this MCMC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_prior_options.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (4) MCMC seetings 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_mcmc_options.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (5) Simulation Initialization, no burnin
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_init.sim.R')


#############################
# Code required
#############################

setwd(paste(home, "/additional_functions/", sep=""))
source(paste("mcmc.extended.full",function_load_name,sep=""))
source(paste("mcmc.full",function_load_name,sep=""))
source(paste("extended_full_log_posterior",function_load_name,sep=""))
source(paste("full_log_posterior",function_load_name,sep=""))
source("logit.R")
source("inv.logit.R")



if (spectral_model=="full"){
  fixedk_mcmc <- mcmc.full
  source("spectral_post_full.R")
} else {
  fixedk_mcmc <- mcmc.extended.full
  source("spectral_post_extended_full.R")
}

setwd(home)



#############################
# Main MCMC run
#############################

online_ordering <- "reference"
no_guess <- k_curr
mu_guess <- mu_curr[,order(w[-mix_num],decreasing=TRUE)]
adapt_end <- 0  # Iteration number on which to stop adaptive MCMC location updates

paras <- vector('list', main_mcmc_iters)
store_log_post <- vector('list', main_mcmc_iters)
store_log_like <- vector('list', main_mcmc_iters)
store_alloc <- allocate_curr*0   # Note, keep this only for no burnin version
store_confusion <- vector('list', main_mcmc_iters)



for (tt in 1:main_mcmc_iters){
  
  # Print iteration at intervals
  if (tt/print_interval == round(tt/print_interval)){
    print(paste("Iteration number: ",tt,sep=""))
    print(paste("Number of sources: ",k_curr,sep=""))
  }
  
  main_run <- fixedk_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas,ewt,k_curr,mix_num)

  # extract allocation and add it to the overall allocation to get allocation rates
  store_alloc = store_alloc + main_run[[2]]

  
  # Extract final parameters and allocation matrix
  new_paras <- main_run[[1]]
  k_curr <- new_paras[1]
  mix_num <- k_curr+1
  mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
  w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]
  eparas <- matrix(new_paras[(3*k_curr+3):(7*k_curr+2)],nrow=k_curr)  
  ewt <- new_paras[(7*k_curr+3):(8*k_curr+2)]
  allocate_curr <- main_run[[2]]
  
  predicted = factor((allocate_curr %*% c(1, 2, 3))[,1], levels = 1:3)
  confusion = table(predicted, truth = dat[,5])

  paras[[tt]] <- new_paras
  store_log_post[[tt]] <- main_run[[3]][1]
  store_log_like[[tt]] <- main_run[[3]][2]
  store_confusion[[tt]] <- confusion
}

main_run <-  list(paras,allocate_curr,store_log_post,store_log_like, store_confusion, alloc_rate = store_alloc/main_mcmc_iters)

# Save output
save(list=ls(),file=paste(results,out.file.name,sep=""))

# Show time taken
end <- proc.time()
print(end-start)

