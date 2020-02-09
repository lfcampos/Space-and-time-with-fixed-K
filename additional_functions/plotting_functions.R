# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Time v Energy distribution funciton
# input time vector, evergy vector for photons
# an object out of kde2d for plotting heatmap
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

plot.TimeEnergy = function(time, energy, k = NULL){
	# define breaks
	breaks = seq(min(time), max(time), length.out = 15)
	# calculate rates
	lambda = dirichlet.mle(time, breakpoints = breaks)
	# plot configuration
	par(mar = c(4, 5, 1, 1), cex.lab=1.5)
	layout(matrix(c(1, 1, 2, 1, 1, 2, 3, 3, 4), byrow = FALSE, nrow = 3))
	# heatplot
	image(k, col =r, ylab = 'log10(Energy)')
	abline(v = breaks, col = 'lightgrey')
	text((breaks[-1] + breaks[-15])/2, log10(rep(min(energy), 14)) + 0.1, 1:14)
	# arrival rate distribution
	draw.bb(breakpoints = breaks, lambda = lambda$N*lambda$alpha, lwd = 2, col = cols[1], xlim = c(min(time) + 900, max(time) - 900), ylab = 'Counts'); abline(v = breaks, col = 'lightgrey')
	# Overall energy distribution
	barplot(hist(log10(energy), plot = FALSE, breaks = 50)$counts, main = '', ylab = 'log10(Energy)', xlab = 'Counts', horiz=T)
	plot(1, type="n", axes=F, xlab="", ylab="")
	
}


plot.TimeEnergy_seperate = function(time, energy, breaks = NULL, k = NULL){
	# calculate rates
	lambda = dirichlet.mle(time, breakpoints = breaks)
	# plot configuration
	par(mar = c(4, 5, 1, 1), cex.lab=1.5)
	layout(matrix(c(1, 1, 2, 1, 1, 2, 3, 3, 4), byrow = FALSE, nrow = 3))
	# heatplot
	for(i in 1:length(k)){
		image(k[[i]], add = (i > 1), xlim = range(breaks), ylim = range(log10(energy)), col =r, ylab = 'log10(Energy)')
	}
	abline(v = breaks, col = 'lightgrey')
	text((breaks[-1] + breaks[-15])/2, log10(rep(min(energy), 14)) + 0.1, 1:14)
	# arrival rate distribution
	draw.bb(breakpoints = breaks, lambda = lambda$N*lambda$alpha, lwd = 2, col = cols[1], xlim = c(min(time) + 900, max(time) - 900), ylab = 'Counts'); abline(v = breaks, col = 'lightgrey')
	# Overall energy distribution
	barplot(hist(log10(energy), plot = FALSE, breaks = 50)$counts, main = '', ylab = 'log10(Energy)', xlab = 'Counts', horiz=T)
	plot(1, type="n", axes=F, xlab="", ylab="")
	
}


