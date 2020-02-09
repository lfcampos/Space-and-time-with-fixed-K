# Extended BASCS algorithm: MCMC for separating overlapping sources
# David Jones, Harvard University
# Luis Campos, Harvard University

####################################################################################
# FIXED K
# NO RJMCMC VERSION - ignore RJMCMC related notes
####################################################################################

####################################################################################
# MARGINAL ENERGY DISTRIBUTION
# Energy Distribution does not vary with time bins
####################################################################################



args=(commandArgs(TRUE))
eval(parse(text = args))

file.name = gsub('./sim_data/', '', file.name)


start <- proc.time()



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

out.file.name = gsub('simulation', 'analysis_eBASCS', file.name)
load(paste(home, data_dir, file.name, sep = ''))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (1) Split data into component pieces and calculate image parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_clean_data.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (2) Set model options, number of breaks, psf, etc.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_model_options.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (3) Prior parameter settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_prior_options.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (4) MCMC seetings 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_mcmc_options.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# (5) Simulation Initialization, no burnin
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source('./additional_functions/_init.time.sim.R')




#############################
# Code required
#############################

setwd(paste(home, "/additional_functions/", sep=""))
source(paste("mcmc.full",function_load_name,sep=""))
source(paste("mcmc.spatial.time.energy.full.marginal",function_load_name,sep=""))
source(paste("mcmc.extended.full",function_load_name,sep=""))
source(paste("mcmc.spatial.time.energy.extended.full.marginal",function_load_name,sep=""))
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
  fixedk_eBASCS_mcmc <- mcmc.spatial.time.energy.extended.full.marginal
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
store_alloc <- allocate_curr*0



for (tt in 1:main_mcmc_iters){
  
  # Print iteration at intervals
  if (tt/print_interval == round(tt/print_interval)){
    print(paste("Iteration number: ",tt,sep=""))
  }
  
  main_run_time <- fixedk_eBASCS_mcmc(mcmc_runs,online_ordering,tt,w,allocate_curr,mu_curr,eparas_all,ewt_all,k_curr,mix_num,bk,num_time_breaks, lambda, time_bin)


  # extract allocation and add it to the overall allocation to get allocation rates
  store_alloc = store_alloc + main_run_time[[2]]

  # Extract parameters and allocation matrix after burnin
  new_paras <- main_run_time[[1]]
  k_curr <- new_paras[1]
  mix_num <- k_curr+1
  mu_curr <- matrix(new_paras[2:(2*k_curr+1)],2,k_curr)
  w <- new_paras[(2*(k_curr+1)):(3*k_curr+2)]

  eparas_all <- main_run_time[[3]] 
  ewt_all <- main_run_time[[4]] 
  eparas_all <- main_run_time[[3]] 
  ewt_all <- main_run_time[[4]] 
  lambda <- main_run_time[[5]] 
  time_bin <- main_run_time[[6]]
  bk <- main_run_time[[7]]
  num_time_breaks <- main_run_time[[8]]

  allocate_curr <- main_run_time[[2]]


  predicted = factor((allocate_curr %*% c(1, 2, 3))[,1], levels = 1:3)
  confusion = table(predicted, truth = dat[,5])

  paras[[tt]] <- list(par = new_paras, eparas_all = eparas_all, ewt_all = ewt_all,lambda = lambda, confusion = confusion)

}
  
main_run_time <- list(paras, allocate_curr, 
  alloc_rate = store_alloc/main_mcmc_iters)

# Save output
save(list=ls(),file=paste(results,out.file.name,sep=""))

# Show time taken
end <- proc.time()
print(end-start)

