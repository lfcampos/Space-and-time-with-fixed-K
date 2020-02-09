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
