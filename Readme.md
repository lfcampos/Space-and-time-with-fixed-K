# Extended BASCS algorithm: MCMC for separating overlapping sources

Luis F. Campos and David E. Jones


**This is an extension of the model described in**: 
``Disentangling Overlapping Astronomical Sources using Spatial and Spectral Information" Jones, Kashyap, van Dyk (ApJ 2015). 

**Follow the writeup for model details in Overleaf**


# Description

The algorithm uses the spatial, temporal and spectral information of a fixed number of sources to find the joint posterior distribution of their spatial, temporal and spectral parameters.

# File Descriptions

- `simulate_data_two_sources.R`: contains code for simulating data from the model for two sources. You can change the distance between them, relative intensities and background level. On the energy side, you we've pre-specified four settings of interest
	1. `equalE`: Equal Energy distributions Gamma(4, 4/1000) across time and sources
	2. `varyE_bySource`: Varying Energy distributions Gamma(4, 4/1000) and Gamma(3, 3/600) for sources, same dist'n across time
	3. `varyE_byTime_Source`: Vary Energy dist'ns across source, bright typically higher than the dim, but vary both across time as well Gamma(~5, ~5/~1000) and Gamma(~3, ~3/~600)
	4. `varyE_byTime_Source_same_Marginal`: Vary Energy dist'n across time for bright source, keep same for dim source -- but marginally will be the ~same. 
The file can be run using the command shell (with optional settings) or manually. 


- Three files to run three different models
	1. `run_BASCS_fixedk.R`: runs BASCS with a fixed number of sources (k), so BASCS without reversible jump. 
	2. `run_eBASCS_full_fixedk.R`: runs eBASCS described in the paper. The energy distributions vary across sources and across time.
	3. `run_eBASCS_full_marginal_fixedk.R`: runs a simplified version of eBASCS where we do not allow the energy distribution to vary across time, but it may vary across sources. 

These files source a number of setting files. We did this to simplify the coding and ensure runs across models have the same settings. These are

- `_clean_data.R`: This formats the data to a form appropriate for the `mcmc` functions
- `_init.sim.R`: Initializes parameters for BASCS mcmc 
- `_init.time.sim.R`: Initializes parameters for eBASCS mcmc 
- `_mcmc_options.R`: Sets parameters for mcmc including the proposal distribution parameters, chain length and burn-in period
- `_model_options.R`: sets model parameter options including spectral model (full or extended), time spacing model (multinomial Dirichlet, equally spaced), psf model. 
- `_prior_options.R`: sets parameters for prior distributions on energy, spatial, etc. 
