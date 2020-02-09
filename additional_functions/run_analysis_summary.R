
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plotting Function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


run_analysis_summary = function(dat, BASCS_run, params = NULL, inc.time = TRUE, make.fig = TRUE, to.pdf = TRUE){

	samples = length(BASCS_run)
	# extract data information
	spatial <- dat[,1:2]
	arrival_time <- dat[,4]
	num <- 2

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Source positions
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	mu_star = t(sapply(BASCS_run[[1]], function(x){x[2:5]}))

	if(make.fig){
		par(mfrow = c(1, 1))
		plot(spatial,main="Two Overlapping Sources", pch=19,cex=0.1)

		plot(NA,main="Posterior: Source centers", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')
		points(mu_star[,1], mu_star[,2], pch = 19,cex=0.2, col = src.col[1])
		points(mu_star[,3], mu_star[,4], pch = 19,cex=0.2, col = src.col[2])

		points(params$positions, pch = 4, cex = 5, lwd = 2)
		legend('topright', src.lab[1:2], col = src.col[1:2], pch = 19)

		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))

	}
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Light Curve posteriors
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


	if(inc.time){
		alloc_source = lapply(BASCS_run, function(x)  x[[2]]%*%1:3)
		bk_source = lapply(BASCS_run, function(x) params$breakpoints)
		lambda_source = vector('list', samples)

		for(i in 1:samples){
		  lambda_source[[i]] = list()
		  for(s in 1:(num+1)){
		    t0 = arrival_time[which(alloc_source[[i]] == s)]
		    if(length(t0) == 0){
		      lambda_source[[i]][[s]] = multinomial.counts(t0, breakpoints = bk_source[[i]][[s]])/length(t0)
		    }else{
		      lambda_source[[i]][[s]] = dirichlet.mle(t0, breakpoints = bk_source[[i]][[s]])$alpha
		    }
		  }
		}


		post_summaries_lambda = list()

		for(s in 1:(num+1)){
		  lambda_source0 = lapply(lambda_source, function(x) x[[s]])
		  lambda_mat = do.call('rbind', lambda_source0)

		  l_mean = apply(lambda_mat, 2, mean, na.rm = TRUE)
		  # l_sd = apply(lambda_mat, 2, sd)
		  l_lower = apply(lambda_mat, 2, quantile, 0.025, na.rm = TRUE)
		  l_upper = apply(lambda_mat, 2, quantile, 0.975, na.rm = TRUE)
		  post_summaries_lambda[[s]] = list(mean = l_mean, lower = l_lower, upper = l_upper)

		  if(make.fig){
		  ylim = range(c(unlist(lambda_source0), params$lambdas[[s]]), na.rm = TRUE)
		  draw.bb(bk_source[[1]][[s]], lambda_source0[[1]], new.plot = TRUE, col = '#d3d3d3', main = paste('Posterior: Light Curves for', src.lab[s]), ylim = ylim)

		  for(r in sample(2:samples, 500)){
		    draw.bb(bk_source[[r]][[s]], lambda_source0[[r]], new.plot = FALSE, col = '#d3d3d3', ylim = ylim)
		  }
		  # draw mean and credible intervals
		  draw.bb(bk_source[[1]][[s]], l_mean, new.plot = FALSE, col = 2, lwd = 1)
		  draw.bb(bk_source[[1]][[s]], l_lower, new.plot = FALSE, col = 2, lty = 3)
		  draw.bb(bk_source[[1]][[s]], l_upper, new.plot = FALSE, col = 2, lty = 3)
		  # draw target
		  draw.bb(params$breakpoints[[s]], params$lambda.dat[[s]], new.plot = FALSE, col = 'blue', lwd = 2)
		  }
		  if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
		}
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Separation and Source Weights
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	dist = sqrt((mu_star[,1] - mu_star[,3])^2 + (mu_star[,2] - mu_star[,4])^2)
	if(make.fig){
		hist(dist, xlab = 'Separation', main = 'Posterior: Separation')
		abline(v = params$separation, col = 'green', lwd = 2)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}


	w_star = t(sapply(BASCS_run[[1]], function(x){x[6:8]}))

	pi.data = (params$obs_nums)/sum(params$obs_nums)

	if(make.fig){
		par(mfrow = c(1, 3), oma=c(0,0,2,0)); 
		for(i in 1:3){
			hist(w_star[,i], main = src.lab[i], xlim = range(c(w_star[,i], pi.data[i])), xlab = expression(w[eval(s)])); abline(v = pi.data[i], col = 'green', lwd = 2)
		}
		title("Posterior: Source Weight", outer=TRUE)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# allocation summaries
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	if(FALSE){
		alloc_source_star = do.call('cbind', alloc_source)
		params.source = dat[,5]
		alloc_corect_star = apply(alloc_source_star, 2, function(x) x == params.source)
		alloc_pct_corect_star = apply(alloc_corect_star, 1, mean)
		par(mfrow = c(1, 3), oma=c(0,0,2,0)); 
		for(i in 1:3){
			hist(alloc_pct_corect_star[params.source == i], main = src.lab[i])
		}
		title("Posterior: Percent Correct Allocations", outer=TRUE)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# calculate means and intervals
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	post_mean_dist = mean(dist)
	post_CI_dist = quantile(dist, c(0.025, 0.975))

	post_mean_w = colMeans(w_star)
	post_CI_w = apply(w_star, 2, quantile, c(0.025, 0.975))
	post_mean_mu = colMeans(mu_star)
	post_CI_mu = apply(mu_star, 2, quantile, c(0.025, 0.975))

	# post_mean_alloc_pct_corect = mean(alloc_pct_corect_star)
	# post_CI_alloc_pct_corect = quantile(alloc_pct_corect_star, c(0.025, 0.975))

	# post_mean_alloc_pct_corect_bySource = tapply(alloc_pct_corect_star, params.source, mean)
	# post_mean_quantile_pct_corect_bySource = tapply(alloc_pct_corect_star, params.source, quantile, c(0.025, 0.975))
	# post_mean_quantile_pct_corect_bySource = tapply(alloc_pct_corect_star, params.source, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))


	post_summaries = list(post_mean_dist = post_mean_dist, 
		post_CI_dist = post_CI_dist, 
		post_mean_w = post_mean_w, 
		post_CI_w = post_CI_w, 
		post_mean_mu = post_mean_mu, 
		post_CI_mu = post_CI_mu)
		# post_mean_alloc_pct_corect = post_mean_alloc_pct_corect, 
		# post_CI_alloc_pct_corect = post_CI_alloc_pct_corect, 
		# post_mean_alloc_pct_corect_bySource = post_mean_alloc_pct_corect_bySource, 
		# post_mean_quantile_pct_corect_bySource = post_mean_quantile_pct_corect_bySource)

	params_keep = list(sep = params$separation, light_curve = diff(params$lambdas[[1]])[1], energy = params$energy_paras, pi_N = params$pi_N, positions = params$positions)
	if(inc.time){
		post_summaries = c(post_summaries, post_summaries_lambda)
	}
	c(post_summaries, params = params_keep)
}













run_analysis_summary_data = function(dat, BASCS_run, inc.time = TRUE, inc.energy = TRUE, make.fig = TRUE, to.pdf = TRUE){

	samples = length(BASCS_run)
	# extract data information
	spatial <- dat[,1:2]
	arrival_time <- dat[,4]
	num <- 2

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Source positions
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	mu_star = t(sapply(BASCS_run[[1]], function(x){x[2:5]}))

	if(make.fig){
		par(mfrow = c(1, 1))
		plot(spatial,main="Posterior Means of Two Overlapping Sources", pch=19,cex=0.2)

		# plot(NA,main="Posterior: Source centers", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')
		points(mu_star[,1], mu_star[,2], pch = 19,cex=0.2, col = src.col[1])
		points(mu_star[,3], mu_star[,4], pch = 19,cex=0.2, col = src.col[2])

		# points(params$positions, pch = 4, cex = 5, lwd = 2)
		legend('topright', src.lab[1:2], col = src.col[1:2], pch = 19)

		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))

	}
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Light Curve posteriors
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


	plot(spatial,main="Final Allocation to Two Overlapping Sources", pch=19,cex=0.2, col = BASCS_run[[2]]%*%c(2, 3, 1))

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Separation and Source Weights
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	dist = sqrt((mu_star[,1] - mu_star[,3])^2 + (mu_star[,2] - mu_star[,4])^2)
	if(make.fig){
		hist(dist, xlab = 'Separation', main = 'Posterior: Separation')
		abline(v = mean(dist), col = 2, lwd = 2)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}


	w_star = t(sapply(BASCS_run[[1]], function(x){x[6:8]}))

	# pi.data = (params$obs_nums)/sum(params$obs_nums)

	if(make.fig){
		par(mfrow = c(1, 3), oma=c(0,0,2,0)); 
		for(i in 1:3){
			hist(w_star[,i], main = src.lab[i], xlim = range(c(w_star[,i])), xlab = expression(w[eval(s)])); abline(v = mean(w_star[,i]), col = 2, lwd = 2)
		}
		title("Posterior: Relative Brightness", outer=TRUE)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))

		par(mfrow = c(1, 1)); 
		plot(density(w_star[,3]), xlim = c(0, 1), col = src.col[3], main = 'Posterior: Relative Brightness', lwd = 2, xlab = 'Brightness')
		lines(density(w_star[,2]), xlim = c(0, 1), col = src.col[2], lwd = 2)
		lines(density(w_star[,1]), xlim = c(0, 1), col = src.col[1], lwd = 2)

		abline(v = apply(w_star, 2, mean), lwd = 1, lty = 2, col = 2)
		legend('topright', legend = src.lab, col = src.col, lwd = 2)

		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))


	}



	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Energy Distribution
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	alloc_source = BASCS_run[[2]]%*%1:3
	
	energy_params = lapply(BASCS_run[[1]], function(x){
		eparas <- matrix(x[(3*k_curr+3):(7*k_curr+2)],nrow=k_curr)  
  		ewt <- x[(7*k_curr+3):(8*k_curr+2)]
  		list(eparas = eparas, ewt = ewt)
  	})

	xlim = c(0, 10000)
	ylim = c(0, 0.0015)
	cols = c(rgb(228/255, 26/255, 28/255, alpha = 0.2), rgb(22/255, 49/255, 72/255, alpha = 0.2))
	par(mfrow = c(1, 1), oma=c(0,0,2,0)); 
	e.grid = seq(0, max(energy), 20)
	for(i in 1:2){	
		hist(dat$pi[alloc_source == i], breaks = 100, prob = TRUE, xlab = 'Energy', main = src.lab[i], xlim = xlim, ylim = ylim)
		tmp = lapply(energy_params[sample(1:length(energy_params), 50)], function(x){
			attach(x)
			lines(e.grid, (ewt[i]*dgamma(e.grid,eparas[i,3],eparas[i,3]/eparas[i,1])+(1-ewt[i])*dgamma(e.grid,eparas[i,4],eparas[i,4]/eparas[i,2])), col = cols[i])
			detach(x)
			return(NULL)
		})
	title("Posterior: Energy Distributions", outer=TRUE)
	}

	return(list(energy_params = energy_params, energy = energy, alloc_source = alloc_source, mu_star = mu_star, w_star = w_star))
}





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# run_analysis_extended_summary_data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #



run_analysis_extended_summary_data = function(dat, BASCS_run, inc.time = TRUE, inc.energy = TRUE, make.fig = TRUE, to.pdf = TRUE){

	samples = length(BASCS_run[[1]])
	# extract data information
	spatial <- dat[,1:2]
	arrival_time <- dat[,4]
	energy <- dat[,3]

	num <- 2

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Source positions
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	mu_star = t(sapply(BASCS_run[[1]], function(x){x$par[2:5]}))

	if(make.fig){
		par(mfrow = c(1, 1))
		plot(spatial,main="Posterior Means of Two Overlapping Sources", pch=19,cex=0.2)

		# plot(NA,main="Posterior: Source centers", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')
		points(mu_star[,1], mu_star[,2], pch = 19,cex=0.2, col = src.col[1])
		points(mu_star[,3], mu_star[,4], pch = 19,cex=0.2, col = src.col[2])

		# points(params$positions, pch = 4, cex = 5, lwd = 2)
		legend('topright', src.lab[1:2], col = src.col[1:2], pch = 19)

		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))

	}
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Light Curve posteriors
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


	plot(spatial,main="Final Allocation to Two Overlapping Sources", pch=19,cex=0.2, col = BASCS_run[[2]]%*%c(2, 3, 1))

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Separation and Source Weights
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	dist = sqrt((mu_star[,1] - mu_star[,3])^2 + (mu_star[,2] - mu_star[,4])^2)
	if(make.fig){
		hist(dist, xlab = 'Separation', main = 'Posterior: Separation')
		abline(v = mean(dist), col = 2, lwd = 2)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}


	w_star = t(sapply(BASCS_run[[1]], function(x){x$par[6:8]}))

	# pi.data = (params$obs_nums)/sum(params$obs_nums)

	if(make.fig){
		par(mfrow = c(1, 3), oma=c(0,0,2,0)); 
		for(i in 1:3){
			hist(w_star[,i], main = src.lab[i], xlim = range(c(w_star[,i])), xlab = expression(w[eval(s)])); abline(v = mean(w_star[,i]), col = 2, lwd = 2)
		}
		title("Posterior: Relative Brightness", outer=TRUE)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))

		par(mfrow = c(1, 1)); 
		plot(density(w_star[,3]), xlim = c(0, 1), col = src.col[3], main = 'Posterior: Relative Brightness', lwd = 2, xlab = 'Brightness')
		lines(density(w_star[,2]), xlim = c(0, 1), col = src.col[2], lwd = 2)
		lines(density(w_star[,1]), xlim = c(0, 1), col = src.col[1], lwd = 2)

		abline(v = apply(w_star, 2, mean), lwd = 1, lty = 2, col = 2)
		legend('topright', legend = src.lab, col = src.col, lwd = 2)

		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))


	}




	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Energy Distributions
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	alloc_source = BASCS_run[[2]] %*% c(1, 2, 3)
	energy_params = lapply(BASCS_run[[1]], function(x){
		eparas <- x$eparas_all 
  		ewt <- x$ewt_all
  		list(eparas = eparas, ewt = ewt)
  	})

	xlim = c(0, 6000)
	ylim = c(0, 0.0015)
	cols = c(rgb(228/255, 26/255, 28/255, alpha = 0.2), rgb(22/255, 49/255, 72/255, alpha = 0.2))
	colfunc <- colorRampPalette(c("springgreen","royalblue"), alpha = 0.2)
	cols = colfunc(14)

	par(mfrow = c(1, 1), oma=c(0,0,2,0)); 
	e.grid = seq(0, max(energy), 20)
	for(i in 1:2){	
		if(sum(alloc_source == i)>0){
		hist(dat$pi[alloc_source == i], breaks = 100, prob = TRUE, xlab = 'Energy', main = src.lab[i], xlim = xlim, ylim = ylim, border = 0)
		tmp = lapply(energy_params[sample(1:length(energy_params), 50)], function(x){
			for(k in 1:length(x[[1]])){
				ewt = x$ewt[[k]]
				eparas = x$eparas[[k]]
				lines(e.grid, (ewt[i]*dgamma(e.grid,eparas[i,3],eparas[i,3]/eparas[i,1])+(1-ewt[i])*dgamma(e.grid,eparas[i,4],eparas[i,4]/eparas[i,2])), col = cols[k])
				
			}
			return(NULL)
		})
		}
	title("Posterior: Energy Distributions", outer=TRUE)
	}

	return(list(energy_params = energy_params, energy = energy, alloc_source = alloc_source, mu_star = mu_star, w_star = w_star))
}







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Summarize recent iteration of ouput 
# May 2018
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


run_analysis_summary_eBASCS = function(dat, BASCS_run, params = NULL, make.fig = TRUE, to.pdf = TRUE, inc.time = FALSE){

	samples = length(BASCS_run)
	# extract data information
	spatial <- dat[,1:2]
	arrival_time <- dat[,4]
	num <- params$num

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Source positions
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	mu_star = t(sapply(BASCS_run[[1]], function(x){x$par[2:5]}))

	if(make.fig){
		par(mfrow = c(1, 1))
		plot(NA,main="Posterior: Source centers", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')
		points(mu_star[,1], mu_star[,2], pch = 19,cex=0.2, col = src.col[1])
		points(mu_star[,3], mu_star[,4], pch = 19,cex=0.2, col = src.col[2])

		points(params$positions, pch = 4, cex = 5, lwd = 2)
		legend('topright', src.lab[1:2], col = src.col[1:2], pch = 19)

		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))

	}
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Light Curve posteriors
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


	if(inc.time){
		lambda_source = list()
		for(i in 1:(num+1)){
			lambda_source[[i]] = lapply(BASCS_run[[1]], function(x)  x$lambda[[i]])
		}
		bk_source = params$breakpoints

		post_summaries_lambda = list()

		for(s in 1:(num+1)){
		  lambda_source0 = lambda_source[[s]]
		  lambda_mat = do.call('rbind', lambda_source0)

		  samples = length(lambda_source0)

		  l_mean = apply(lambda_mat, 2, mean, na.rm = TRUE)
		  # l_sd = apply(lambda_mat, 2, sd)
		  l_lower = apply(lambda_mat, 2, quantile, 0.025, na.rm = TRUE)
		  l_upper = apply(lambda_mat, 2, quantile, 0.975, na.rm = TRUE)
		  post_summaries_lambda[[s]] = list(mean = l_mean, lower = l_lower, upper = l_upper)

		  if(make.fig){
		  ylim = range(c(unlist(lambda_source0), params$lambdas[[s]]), na.rm = TRUE)
		  draw.bb(bk_source[[s]], lambda_source0[[1]], new.plot = TRUE, col = '#d3d3d3', main = paste('Posterior: Light Curves for', src.lab[s]), ylim = ylim)

		  for(r in sample(2:samples, 500)){
		    draw.bb(bk_source[[s]], lambda_source0[[r]], new.plot = FALSE, col = '#d3d3d3', ylim = ylim)
		  }
		  # draw mean and credible intervals
		  draw.bb(bk_source[[s]], l_mean, new.plot = FALSE, col = 2, lwd = 1)
		  draw.bb(bk_source[[s]], l_lower, new.plot = FALSE, col = 2, lty = 3)
		  draw.bb(bk_source[[s]], l_upper, new.plot = FALSE, col = 2, lty = 3)
		  # draw target
		  draw.bb(params$breakpoints[[s]], params$lambda.dat[[s]], new.plot = FALSE, col = 'blue', lwd = 2)
		  }
		  if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
		}
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# Separation and Source Weights
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	dist = sqrt((mu_star[,1] - mu_star[,3])^2 + (mu_star[,2] - mu_star[,4])^2)
	if(make.fig){
		hist(dist, xlab = 'Separation', main = 'Posterior: Separation')
		abline(v = params$separation, col = 'green', lwd = 2)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}


	w_star = t(sapply(BASCS_run[[1]], function(x){x$par[6:8]}))

	pi.data = (params$obs_nums)/sum(params$obs_nums)

	if(make.fig){
		par(mfrow = c(1, 3), oma=c(0,0,2,0)); 
		for(i in 1:3){
			hist(w_star[,i], main = src.lab[i], xlim = range(c(w_star[,i], pi.data[i])), xlab = expression(w[eval(s)])); abline(v = pi.data[i], col = 'green', lwd = 2)
		}
		title("Posterior: Source Weight", outer=TRUE)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# allocation summaries
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


	if(FALSE){
	lapply(BASCS_run[[1]], function(x) {
		x$confusion
		
	})
		alloc_source_star = do.call('cbind', alloc_source)
		params.source = dat[,5]
		alloc_corect_star = apply(alloc_source_star, 2, function(x) x == params.source)
		alloc_pct_corect_star = apply(alloc_corect_star, 1, mean)
		par(mfrow = c(1, 3), oma=c(0,0,2,0)); 
		for(i in 1:3){
			hist(alloc_pct_corect_star[params.source == i], main = src.lab[i])
		}
		title("Posterior: Percent Correct Allocations", outer=TRUE)
		if(!to.pdf) invisible(readline(prompt="Press [enter] to continue"))
	}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
	# calculate means and intervals
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

	post_mean_dist = mean(dist)
	post_CI_dist = quantile(dist, c(0.025, 0.975))

	post_mean_w = colMeans(w_star)
	post_CI_w = apply(w_star, 2, quantile, c(0.025, 0.975))
	post_mean_mu = colMeans(mu_star)
	post_CI_mu = apply(mu_star, 2, quantile, c(0.025, 0.975))

	# post_mean_alloc_pct_corect = mean(alloc_pct_corect_star)
	# post_CI_alloc_pct_corect = quantile(alloc_pct_corect_star, c(0.025, 0.975))

	# post_mean_alloc_pct_corect_bySource = tapply(alloc_pct_corect_star, params.source, mean)
	# post_mean_quantile_pct_corect_bySource = tapply(alloc_pct_corect_star, params.source, quantile, c(0.025, 0.975))
	# post_mean_quantile_pct_corect_bySource = tapply(alloc_pct_corect_star, params.source, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))


	post_summaries = list(post_mean_dist = post_mean_dist, 
		post_CI_dist = post_CI_dist, 
		post_mean_w = post_mean_w, 
		post_CI_w = post_CI_w, 
		post_mean_mu = post_mean_mu, 
		post_CI_mu = post_CI_mu)
		# post_mean_alloc_pct_corect = post_mean_alloc_pct_corect, 
		# post_CI_alloc_pct_corect = post_CI_alloc_pct_corect, 
		# post_mean_alloc_pct_corect_bySource = post_mean_alloc_pct_corect_bySource, 
		# post_mean_quantile_pct_corect_bySource = post_mean_quantile_pct_corect_bySource)

	params_keep = list(sep = params$separation, light_curve = diff(params$lambdas[[1]])[1], energy = params$energy_paras, pi_N = params$pi_N, positions = params$positions)

	if(inc.time){
		post_summaries = c(post_summaries, post_summaries_lambda)
	}
	c(post_summaries, params = params_keep)
}
