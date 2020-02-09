
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Rewriting the Vanderplas code in R
# Based on algorithm outlined in http://adsabs.harvard.edu/abs/2012arXiv12
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
bayesian_blocks = function(t, p0 = 0.01){
	t = sort(t)
	N = length(t)
	
	# formula (21) for penalty from prior: see Studies in Astronomical Time Series Analysis. Scragle - http://arxiv.org/abs/1207.5578
	prior = 4 - 73.53*p0*N^{-0.478}

	# create length-(N + 1) array of cell edges
	edges = c(t[1], 0.5*(t[2:N] + t[1:(N-1)]), t[N])
	block_length = t[N] - edges

	# arrays needed for the iteration
	nn_vec = rep(1, N)
	best = rep(0, N)
	last = rep(0, N)

	#-----------------------------------------------------------------
	# Start with first data cell; add one cell at each iteration
	#-----------------------------------------------------------------

	for(K in 1:N){
		# Compute the width and count of the final bin for all possible
		# locations of the K^th changepoint
		width = block_length[1:K] - block_length[K + 1]
		count_vec = rev(cumsum(rev(nn_vec[1:K])))
 
		# evaluate fitness function for these possibilities
		fit_vec = count_vec*(log(count_vec) - log(width))
		fit_vec = fit_vec - prior #PRIOR
		fit_vec[2:K] = fit_vec[2:K] + best[1:(K-1)]

		i_max = which.max(fit_vec)
		last[K] = i_max
		best[K] = fit_vec[i_max]
	}


	change_points = rep(0, N)
	i_cp = N + 1
	ind = N + 1

	while(TRUE){
		i_cp = i_cp - 1
		change_points[i_cp] = ind
		if(ind == 1) break()
		ind = last[ind-1]
	}

	edges[change_points[i_cp:N]]
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Light Curve Plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
lc.plot = function(t, breaks = 100, ylim = c(0, 50), main = 'Light Curve', blocks = FALSE, blocks.b = NULL){
	h = hist(t, breaks = breaks, lty = 0, ylim = ylim, main = main)
	points(h$mids, h$counts, pch = 20, type = 'b')
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Fits the Bayesian blocks model with Poisson proccess mean estimate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
fit.bb = function(t, breakpoints = NULL){
	L = diff(breakpoints)
	lambda = table(cut(t, breakpoints))/L
	lambda
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Find MLE of the Dirichlet distribution model given times and breakpoints
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
dirichlet.mle = function(t, breakpoints = NULL){
	N = length(t)
	beta = table(cut(t, breakpoints))/N
	list(N = N, alpha = beta)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Find multinomial counts of times (t) in bins defined by breakpoints
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
multinomial.counts = function(t, breakpoints = NULL){
	table(cut(t, breakpoints))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Draws Bayesian blocks model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
draw.bb = function(breakpoints, lambda, new.plot = TRUE, col = 1, main = '', ylim = NULL, lwd = 1, ylab = NULL, lines = FALSE, lty = 1, xlim = NULL){
	if(is.null(ylim)){ ylim = range(lambda) }
	if(is.null(xlim)){ xlim = range(breakpoints)}
	if(is.null(ylab)){ylab = expression(lambda)}
	n.seg = length(lambda)
	if(new.plot){ plot(NA, xlim = xlim, ylim = ylim, type = 'n', xlab = 't', ylab = ylab, main = main) }
	segments(x0 = breakpoints[1:n.seg], y0 = lambda, x1 = breakpoints[2:(n.seg+1)], col = col, lwd = lwd, lty = lty)
	if(lines){
		points((breakpoints[-1] + breakpoints[-length(breakpoints)])/2, lambda, type = 'b', col = col, lwd = lwd)
	}
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Gets the rates a Poisson Proccesses of with a step funciton 
#	at time t0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

get.rates = function(t0, breaks, lambda){
	# findInterval(t0, breaks)
	lambda[findInterval(t0, breaks, rightmost.closed = TRUE)]
}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Calculates the difference of two step functions as defined by their
#    breakpoints (length n) and 
#    heights (length n-1)
# Returns: new step function, the difference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

lambda_diff = function(b1, b2, lambda1, lambda2){
	# extend last and first step function value to min and max of both breakpoints
	#	essentially aligning the two 
	b1 = sort(b1)
	b2 = sort(b2)
	b1[1] <- b2[1] <- min(c(b1, b2))
	b1[length(b1)] <- b2[length(b2)] <- max(c(b1, b2))

	b = sort(unique(c(b1, b2)))
	cs = (b[-length(b)] + b[-1])/2

	l.diff = get.rates(cs, b1, lambda1) - get.rates(cs, b2, lambda2)

	list(bk = b, lambda.diff = l.diff)
}
