# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : summarize_simulations_inspection.R
# Programmer Name    : Luis Campos
#                     soyluiscampos@gmail.com
#
# Purpose            : We inspect wether the (multiple) chains are 
#                      converging and not getting stuck at local modes
#
# Input              : None
# Output             : None
# Usage              : 
# 
# References         : None
#
#
# Platform           : R
# Version            : v3.3.0
# Date               : 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

source('summarize_functions.R')

names = as.character(read.table(file = './results_chains/filenames.txt')[,1])

ids = do.call('rbind', strsplit(gsub('analysis_eBASCS_varyE_byTime_Source_param_set_|\\.Rda', '', names), '_rep_'))
ids = data.frame(ids)
names(ids) = c('setting', 'chain')


setting = 13

to.load = names[ids$setting == setting]
eBASCS_output = lapply(to.load, function(x){
	load(paste('./results_chains/', x, sep = ''))
	main_run_time
})


posterior = MCMC_matrix_form(eBASCS_output)


par_names = dimnames(posterior)[[3]]

mcmc_trace(posterior, pars = c("mu_1_x", "mu_1_y", "mu_2_x", "mu_2_y"))
mcmc_trace(posterior, pars = c("w_1", "w_2", "w_0"))
mcmc_trace(posterior, pars = par_names[grep('eparas_1', par_names)])
mcmc_trace(posterior, pars = par_names[grep('eparas_2', par_names)])
mcmc_trace(posterior, pars = par_names[grep('lambda_1', par_names)])
mcmc_trace(posterior, pars = par_names[grep('lambda_2', par_names)])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# I've come up with a clean way of visualizing the MCMC output [check]
# I need to 1, understand and show that the chan behavior in the "true"
# blocks case is indeed behaving as it should. Then dig into why this is 
# happening. Run this, same number of chains, given the true blocks
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #




mcmc_dens_overlay(posterior, pars = c("mu_1_x", "mu_1_y", "mu_2_x", "mu_2_y"))
mcmc_dens_overlay(posterior, pars = c("w_1", "w_2", "w_0"))
mcmc_dens_overlay(posterior, pars = par_names[grep('eparas_1', par_names)])
mcmc_dens_overlay(posterior, pars = par_names[grep('eparas_2', par_names)])
mcmc_dens_overlay(posterior, pars = par_names[grep('lambda_1', par_names)])
mcmc_dens_overlay(posterior, pars = par_names[grep('lambda_2', par_names)])

mcmc_areas(posterior, pars = c("mu_1_x", "mu_1_y", "mu_2_x", "mu_2_y"), prob = 0.8, prob_outer = 0.99, point_est = "mean")
mcmc_areas(posterior, pars = c("w_1", "w_2", "w_0"), prob = 0.8, prob_outer = 0.99, point_est = "mean")
mcmc_areas(posterior, pars = par_names[grep('eparas_1', par_names)], prob = 0.8, prob_outer = 0.99, point_est = "mean")
mcmc_areas(posterior, pars = par_names[grep('eparas_2', par_names)], prob = 0.8, prob_outer = 0.99, point_est = "mean")
mcmc_areas(posterior, pars = par_names[grep('lambda_1', par_names)], prob = 0.8, prob_outer = 0.99, point_est = "mean")
mcmc_areas(posterior, pars = par_names[grep('lambda_2', par_names)], prob = 0.8, prob_outer = 0.99, point_est = "mean")



mcmc_pairs(posterior, pars = c("w_1", "w_2", "w_0"))
