adapt_end <- 1000  # Iteration number on which to stop adaptive MCMC location updates
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
burnin_mcmc_iters <- 3000  # Number of RJMCMC iterations for the burnin period
main_mcmc_iters <- 10000 # Number of RJMCMC iterations for the main run

# Thinning:
mcmc_runs <- 10  # Number of steps of MCMC to run between each recording or parameters
print_interval <- 500  # Number of iterations between print of current iteration number 

