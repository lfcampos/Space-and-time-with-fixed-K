# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# R CMD BATCH --no-save --no-restore "--args eng.setting='equalE'" simulate_data_two_sources.R ./simulate_data_two_sources_1.out 
# R CMD BATCH --no-save --no-restore "--args eng.setting='varyE_bySource'" simulate_data_two_sources.R ./simulate_data_two_sources_2.out 
# R CMD BATCH --no-save --no-restore "--args eng.setting='varyE_byTime_Source'" simulate_data_two_sources.R ./simulate_data_two_sources_3.out 
# R CMD BATCH --no-save --no-restore "--args eng.setting='varyE_byTime_Source_same_Marginal'" simulate_data_two_sources.R ./simulate_data_two_sources_4.out 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Settings Descriptions
# separation - distance between two sources: in the range of 0 - 1.5
# rel_intensity - The relative intesity of the bright source to the dim source: > 1
# rel_back - Relative intensity of background with maximal being ~5% of the dim source, this is when rel_back = 1. Kept relatively small.
# datanum - Replicates dataset 10 times to assess variability within settings

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Energy Distribution option (eng.setting): 
# 1. "equalE": Equal Energy distributions Gamma(4, 4/1000) accross time and sources
# 2. "varyE_bySource": Varying Energy distributions Gamma(4, 4/1000) and Gamma(3, 3/600) for sources, same dist'n accross time
# 3. "varyE_byTime_Source": Vary Energy dist'ns accross source, bright typically higher than the dim, but vary both accross time as well Gamma(~5, ~5/~1000) and Gamma(~3, ~3/~600)
# 4. "varyE_byTime_Source_same_Marginal": Vary Energy dist'n accross time for bright source, keep same for dim source -- but marginally will be the ~same. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
args=(commandArgs(TRUE))
eval(parse(text = args))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Manually set the enery settings, can be read as an argument 
eng.setting='equalE'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


library(cubature)


dir_rking = './psf_rand/'
source(paste(dir_rking, "rsking.R", sep = ''))
source(paste(dir_rking, "radius.R", sep = ''))
source(paste(dir_rking, "king_profile_density.R", sep = ''))
source(paste(dir_rking, "king.constant.R", sep = ''))

dir_functions = './additional_functions/'
source(paste(dir_functions, "light_curves.R", sep = ''))


out_fig = './sim_figures/'
out_dat = './sim_data/'


file.name <- NULL

n.bins = 5
bp = seq(0, 60000, length.out = n.bins+1)
breakpoints = list(bp, bp, bp)

lambdas_name = 'linear'
x = seq(-n.bins/2, n.bins/2, length.out = n.bins)
deltas = seq(0, 1/(-x[1]*n.bins + 2), length.out = 2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Parameters: sort right to left then switch (5:1)
# 
# delta: controls time informativeness (0 is flat accross time for bright)
# 		 dim is always flat
# r: Relative Intensity - dim = bright/r
# b: relative background - background = b*dim
# d: distance - bright to dim in standard units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
params_all = expand.grid(datanum = 1, r = c(2), b = c(0.1), d = c(1.5), delta = deltas)[,5:1]

if(eng.setting == "equalE"){
	# Set Energy Parrameters 
	num <- 2 
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][,1] = c(1832, 1832)
		energy_paras[[i]][,2] = c(3.18, 3.18)
	}
}

if(eng.setting == "varyE_bySource"){
	# Set Energy Parrameters 
	num <- 2 
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][,1] = c(1000, 600)
		energy_paras[[i]][,2] = c(4, 3)
	}
}

if(eng.setting == "varyE_byTime_Source"){
	# Set Energy Parrameters 
	num <- 2 
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][1,1] = rnorm(1, 1000, 250)
		energy_paras[[i]][2,1] = rnorm(1, 600, 200)
		energy_paras[[i]][1,2] = rnorm(1, 5, 0.5)
		energy_paras[[i]][2,2] = rnorm(1, 3, 0.5)
	}
}


if(eng.setting == "varyE_byTime_Source_same_Marginal"){
	# Set Energy Parrameters 
	num <- 2 
	add = seq(0, 5, length.out = n.bins)
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][1,1] = 1000*(3 + add[i])/3 # keep gamma/beta constant
		energy_paras[[i]][2,1] = 1832				 # generates data similar to marginal of bright source
		energy_paras[[i]][1,2] = 3 + add[i]
		energy_paras[[i]][2,2] = 3.18                # generates data similar to marginal of bright source
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Parameters for King Profile and Image Boundaries
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# parameters for king Profile
r0 <- 0.6 # 5.128343
slope <- 1.5 # 1.58749
ellip <- 0.00573533
off.angle <- 0

# set boundaries for 
xup <- 5
xlow <- -5
yup <- 5
ylow <- -5
max_back_energy <- 10000 # for uniform background


pb <- txtProgressBar(min = 0, max = nrow(params_all), style = 3)
pbcounter = 1

for(p in 1:nrow(params_all)){

	# Set Counter
	pbcounter = pbcounter + 1
	setTxtProgressBar(pb, pbcounter)
	
	# Extract Parameters
	delta = params_all[p,1]
	separation = params_all[p,2]
	rel_back = params_all[p,3]
	rel_intensity = params_all[p,4]
	datanum = params_all[p,5]

    # Calculate light curve relative intensity
	l1 = sapply(1/n.bins + x*delta, function(x) max(x, 0))
	l1 = l1/sum(l1)
	l2 = rep(1/n.bins, n.bins)
	l3 = rep(1/n.bins, n.bins)
	lambdas = list(l1, l2, l3)


	# Specify number of sources and either enter positions 
	num <- 2 
	positions <- matrix(NA,num,2)
	positions[1,] <- c(-separation/2,0)
	positions[2,] <- c(separation/2,0) 


	# calculate rates for source intensities
	yl <- (yup - ylow)/2
	xl <- (xup - xlow)/2
	pi_N <- matrix(NA,num+1,1)
	pi_N[1] <- 2000
	pi_N[2] <- pi_N[1]/rel_intensity
	pi_N[3] <- rel_back*pi_N[2]*0.513/10.7 


	# arrange light curves data for simulation	
	b1 = breakpoints[[1]]; b2 = breakpoints[[2]]
	lambda1 = lambdas[[1]]; lambda2 = lambdas[[2]]

	# gather light curve information into one list for use later
	light_curves = list('bright' = cbind(min = b1[-length(b1)], 
		max = b1[-1], rate = lambda1), 'dim' = cbind(min = b2[-length(b2)], 
		max = b2[-1], rate = lambda2))

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Generate data
	#	Source intensities Pois(pi_N)
	#	Spatial Data: king distribution with parameters set above
	#	Energy data: keep v1 energy settings (for later)
	#	Arrival Time: Multinomial with 
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	N_sim <- matrix(NA,num+1,1)
	k_bin <- matrix(NA,num,1)
	sall_sim_v1 <- vector("list",num)
	# only keep first energy parameters
	for (i in 1:num){
		N_sim[i] <- rpois(1,pi_N[i])
		king.norm <- 1/king.constant(100)$integral
		## Simulate Location
		idraws <- rsking(N_sim[i],r0,slope)+cbind(rep(positions[i,1],N_sim[i]),rep(positions[i,2],N_sim[i])) 
		## Simulate Arrival Time
		k_bin[i] = nrow(light_curves[[i]])
		# decide which bin of the curve each photon will fall under
		which.bk = sample(1:k_bin[i], N_sim[i], replace = TRUE, prob = light_curves[[i]][,'rate'])
		# once the bin is decided, generate time as a uniform distribution 
		#   within the selected bin
		t = light_curves[[i]][which.bk,'min'] + runif(N_sim[i])*(light_curves[[i]][which.bk,'max'] - light_curves[[i]][which.bk,'min'])


		## Simulate Energy
		alphas_i = sapply(energy_paras[which.bk], function(x) x[i,1])
		gammas_i = sapply(energy_paras[which.bk], function(x) x[i,2])

		## #LUIS HERE: Fix these to match MCMC formulation. 
		eng <- rgamma(N_sim[i], gammas_i, gammas_i/alphas_i)

		sall_sim_v1[[i]] <- cbind(idraws, eng, t)	
		row.names(sall_sim_v1[[i]]) <- NULL
	}


	# Delete observations outside image
	obs_sim_v1 <- vector("list",num)
	for (i in 1:num){
		obs_sim_v1[[i]] <- sall_sim_v1[[i]][sall_sim_v1[[i]][,1] >= xlow & sall_sim_v1[[i]][,1] <= xup & sall_sim_v1[[i]][,2] >= ylow & sall_sim_v1[[i]][,2] <= yup,]
	}
	obs_nums <- matrix(NA,num+1,1)
	for (i in 1:num){
		obs_nums[i] <- length(obs_sim_v1[[i]][,1])
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Simulate Background
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	N_sim[num+1] <- rpois(1,4*xl*yl*pi_N[num+1])
	xb <- runif(N_sim[num+1],xlow,xup)
	yb <- runif(N_sim[num+1],ylow,yup)
	eb_v1 <- runif(N_sim[num+1], 0, max_back_energy)
	# rgamma(N_sim[num+1],energy_paras[1,num+1],energy_paras[2,num+1])
	obs_nums[num+1] <- N_sim[num+1]

	# times data distributed for background
	t.min = min(unlist(breakpoints))
	t.max = max(unlist(breakpoints))
	t0 = runif(N_sim[num+1], t.min, t.max)

	obs_sim_v1_b = cbind(xb, yb, eb_v1, t0)


	lambda.dat = vector('list', num+1)
	for (i in 1:num){
		lambda.dat[[i]] = dirichlet.mle(sall_sim_v1[[i]][,4], breakpoints[[i]])$alpha
	}
	breakpoints[[num+1]] = breakpoints[[1]]
	lambda.dat[[num+1]] = dirichlet.mle(t0, breakpoints[[num+1]])$alpha


	dat = cbind(rbind(do.call('rbind', obs_sim_v1), obs_sim_v1_b), source = rep(1:3, obs_nums))

	params = list(num = num, xup = xup, xlow = xlow, yup = yup, ylow = ylow, 
		separation = separation, rel_intensity = rel_intensity, 
		rel_back = rel_back, positions = positions, pi_N = pi_N, 
		N_sim = N_sim, obs_nums = obs_nums, energy_paras = energy_paras, 
		max_back_energy = max_back_energy, breakpoints = breakpoints, 
		lambdas = lambdas, lambda.dat = lambda.dat)

	fn = paste(out_dat, 'simulation_', eng.setting,'_param_set_', 
		p, '.Rda', sep = '')

	save(dat, params, file = fn)


	file.name = c(file.name, fn)

}

close(pb)


save(params_all, file = paste('params_all.Rda', sep = ''))




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Some plotting functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if(FALSE){
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Some simulations for "varyE_byTime_Source_same_Marginal"
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	plot_eng = function(energy_paras){
		bright_eng = lapply(energy_paras, function(x) rgamma(1000, x[1,2], x[1,2]/x[1,1]))
		dim_eng = lapply(energy_paras, function(x) rgamma(1000, x[2,2], x[2,2]/x[2,1]))
		tmp1 = data.frame(eng = unlist(bright_eng), time = sort(rep(1:n.bins, 1000)))
		tmp2 = data.frame(eng = unlist(dim_eng), time = sort(rep(1:n.bins, 1000)))
		par(mfrow = c(2, 2), mar=c(3.1,4.1,1.5,0.5), oma = c(0,0,3.5,0))
		ylim = range(c(tmp1$eng, tmp2$eng))
		boxplot(eng~time, data = tmp1, ylim = ylim, main = 'Bright Source Energy')
		boxplot(eng~time, data = tmp2, ylim = ylim, col = rgb(1, 0, 0, alpha = 0.2), main = 'Dim Source Energy')
		bk = seq(min(c(tmp1$eng, tmp2$eng)), max(c(tmp1$eng, tmp2$eng)), length.out = 100)
		p1 = hist(tmp1$eng, breaks = bk, plot = FALSE)
		p2 = hist(tmp2$eng, breaks = bk, plot = FALSE)
		hist(tmp1$eng, ylim = c(0, max(c(p1$counts,  p2$counts))), breaks = bk, main = 'Marginals', xlab = 'Energy')
		hist(tmp2$eng, ylim = c(0, max(c(p1$counts,  p2$counts))), breaks = bk, add = TRUE, col = rgb(1, 0, 0, alpha = 0.2))
		legend('topright', c('Bright', 'Dim'), fil = c(rgb(1,1,1), rgb(1,0,0)))

		# # best fit gamma
		# m = glm(eng~1, data = tmp, family = Gamma(link = "identity"))
		# hist(tmp$eng, prob = TRUE, breaks = 100)
		# lines(xtic, dgamma(xtic, gamma.shape(m)[[1]], gamma.shape(m)[[1]]/coefficients(m)))

		heights = lapply(deltas, function(delta){
			sapply(1/n.bins + x*delta, function(x) max(x, 0))
		})
		plot(heights[[1]], ylim = c(0, max(unlist(heights))), type = 'b', col = 2, main = 'Relative Intensities', ylab = expression(lambda), xlab = 'time')
		lapply(heights[-1], function(x){points(x, type = 'b')})
		title(eng.setting, outer = TRUE)
	}



pdf(file = 'Energy Simulation Settings.pdf')

eng.setting = "equalE"
if(eng.setting == "equalE"){
	# Set Energy Parrameters 
	num <- 2 
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][,1] = c(1000, 1000)
		energy_paras[[i]][,2] = c(4, 4)
	}
	plot_eng(energy_paras)
}

eng.setting = "varyE_bySource"
if(eng.setting == "varyE_bySource"){
	# Set Energy Parrameters 
	num <- 2 
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][,1] = c(1000, 600)
		energy_paras[[i]][,2] = c(4, 3)
	}
	plot_eng(energy_paras)
}

eng.setting = "varyE_byTime_Source"
if(eng.setting == "varyE_byTime_Source"){
	# Set Energy Parrameters 
	num <- 2 
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][1,1] = rnorm(1, 1000, 250)
		energy_paras[[i]][2,1] = rnorm(1, 600, 200)
		energy_paras[[i]][1,2] = rnorm(1, 5, 0.5)
		energy_paras[[i]][2,2] = rnorm(1, 3, 0.5)
	}
	plot_eng(energy_paras)

}


eng.setting = "varyE_byTime_Source_same_Marginal"
if(eng.setting == "varyE_byTime_Source_same_Marginal"){
	# Set Energy Parrameters 
	num <- 2 
	add = seq(0, 5, length.out = n.bins)
	energy_paras = list()
	for(i in 1:n.bins){
		energy_paras[[i]] <- matrix(NA,2,num)
		energy_paras[[i]][1,1] = 1000*(3 + add[i])/3 # keep gamma/beta constant
		energy_paras[[i]][2,1] = 1832
		energy_paras[[i]][1,2] = 3 + add[i]
		energy_paras[[i]][2,2] = 3.18
	}
	plot_eng(energy_paras)

}

dev.off()




}