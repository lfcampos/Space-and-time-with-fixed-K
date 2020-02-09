library("bayesplot")
library("ggplot2")
library('Rcpp')
library("rstanarm") 
library('abind')




library('MASS')
library('plotrix')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Clean MCMC output to matrices, maybe clean up for MCMC Visualization form
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


MCMC_matrix_form = function(eBASCS_output){
	# number of chains
	num.chains = length(eBASCS_output)
	# assuming same number of iterations per chain
	num.iter = length(eBASCS_output[[1]][[1]]) 
	#number of parameters per chain

	# gather parameters (list item for each chain)
	par_arrays = lapply(eBASCS_output, function(one_chain) do.call('rbind', lapply(one_chain[[1]], function(x) x$par)))
	for(i in 1:num.chains){
		colnames(par_arrays[[i]]) = c('k', 'mu_1_x', 'mu_1_y', 'mu_2_x', 'mu_2_y', 'w_1', 'w_2', 'w_0')
	}
	epar_arrays = lapply(eBASCS_output, function(one_chain) do.call('rbind', lapply(one_chain[[1]], function(x){ 
			e_list = lapply(x$eparas_all, unlist)
			for(i in 1:length(e_list)){
				x = e_list[[i]]
				names(e_list[[i]]) = paste('eparas_', i, '_', sort(rep(1:(length(x)/2), 2)), rep(1:2, length(x)/2), sep = '')
			}
			unlist(e_list)
		})
	))
	lambda_arrays = lapply(eBASCS_output, function(one_chain) do.call('rbind', lapply(one_chain[[1]], function(x){
			l = x$lambda
			for(i in 1:length(l)){
				l[[i]] = l[[i]][1,]
				names(l[[i]]) = paste('lambda', i, 1:length(l[[i]]), sep = '_')
			}
			unlist(l)
		})
	))

	# organize into array form  [iteration, chain, parameter]
	par_array = array(NA, c(num.iter, num.chains, ncol(par_arrays[[1]])))
	epar_array = array(NA, c(num.iter, num.chains, ncol(epar_arrays[[1]])))
	lambda_array = array(NA, c(num.iter, num.chains, ncol(lambda_arrays[[1]])))

	for(i in 1:num.chains){
		par_array[, i, ] = par_arrays[[i]]
		epar_array[, i, ] = epar_arrays[[i]]
		lambda_array[, i, ] = lambda_arrays[[i]]
	}

	dimnames(par_array) = list(iterations = NULL, chains = paste('chain', 1:num.chains, sep = ':'), parameters = colnames(par_arrays[[1]]))
	dimnames(epar_array) = list(iterations = NULL, chains = paste('chain', 1:num.chains, sep = ':'), parameters = colnames(epar_arrays[[1]]))
	dimnames(lambda_array) = list(iterations = NULL, chains = paste('chain', 1:num.chains, sep = ':'), parameters = colnames(lambda_arrays[[1]]))

	posterior = abind(par_array, epar_array, lambda_array, along = 3)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot some simple Ellipses mean +/- 2*sd and corr rotation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
subset.file = function(which.eng.setting, which.time.setting){
	file.eBASCS[with(file.params, which(eng.setting == which.eng.setting & time.setting == which.time.setting))]
}
library(plotrix)
plot_location_ellipses = function(which.eng.setting, which.time.setting){

	fs = subset.file(which.eng.setting, which.time.setting)

	par(mfrow = c(1, 2), mar=c(3.1,4.1,1.5,0.5), oma = c(0,0,3.5,0))

	plot(NA,main="BASCS", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')

	for(f in fs){
		load(gsub('\\_eBASCS', '', f))
		mu_star_BASCS = t(sapply(main_run[[1]], function(x){x[2:5]}))
		draw.ellipse(mean(mu_star_BASCS[,1]),  mean(mu_star_BASCS[,2]), 2*sd(mu_star_BASCS[,1]),  2*sd(mu_star_BASCS[,2]), cor(mu_star_BASCS[,1], mu_star_BASCS[,2])*180, deg = TRUE, border = src.col[1], lwd = 1.5)
		draw.ellipse(mean(mu_star_BASCS[,3]),  mean(mu_star_BASCS[,4]), 2*sd(mu_star_BASCS[,3]),  2*sd(mu_star_BASCS[,4]), cor(mu_star_BASCS[,3], mu_star_BASCS[,4])*180, deg = TRUE, border = src.col[2], lwd = 1.5)
	}

	plot(NA,main="eBASCS", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')

	for(f in fs){
		load(f)
		mu_star = t(sapply(main_run_time[[1]], function(x){x$par[2:5]}))
		draw.ellipse(mean(mu_star[,1]),  mean(mu_star[,2]), 2*sd(mu_star[,1]),  2*sd(mu_star[,2]), cor(mu_star[,1], mu_star[,2])*180, deg = TRUE, border = src.col[1], lwd = 1.5)
		draw.ellipse(mean(mu_star[,3]),  mean(mu_star[,4]), 2*sd(mu_star[,3]),  2*sd(mu_star[,4]), cor(mu_star[,3], mu_star[,4])*180, deg = TRUE, border = src.col[2], lwd = 1.5)
	}

	title(paste('Energy:', which.eng.setting, '    and     Time:', which.time.setting), outer = TRUE)


}

plot_location_ellipses2 = function(which.eng.setting, which.time.setting, method = 'BASCS', main = ''){

	fs = subset.file(which.eng.setting, which.time.setting)

	plot(NA,main=main, pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = '', ylab = '')

	if(method == "BASCS"){
		for(f in fs){
			load(gsub('\\_eBASCS', '', f))
			mu_star_BASCS = t(sapply(main_run[[1]], function(x){x[2:5]}))
			draw.ellipse(mean(mu_star_BASCS[,1]),  mean(mu_star_BASCS[,2]), 2*sd(mu_star_BASCS[,1]),  2*sd(mu_star_BASCS[,2]), cor(mu_star_BASCS[,1], mu_star_BASCS[,2])*180, deg = TRUE, border = src.col[1], lwd = 1.5)
			draw.ellipse(mean(mu_star_BASCS[,3]),  mean(mu_star_BASCS[,4]), 2*sd(mu_star_BASCS[,3]),  2*sd(mu_star_BASCS[,4]), cor(mu_star_BASCS[,3], mu_star_BASCS[,4])*180, deg = TRUE, border = src.col[2], lwd = 1.5)
		}
	}

	if(method == "eBASCS"){
		for(f in fs){
			load(f)
			mu_star = t(sapply(main_run_time[[1]], function(x){x$par[2:5]}))
			draw.ellipse(mean(mu_star[,1]),  mean(mu_star[,2]), 2*sd(mu_star[,1]),  2*sd(mu_star[,2]), cor(mu_star[,1], mu_star[,2])*180, deg = TRUE, border = src.col[1], lwd = 1.5)
			draw.ellipse(mean(mu_star[,3]),  mean(mu_star[,4]), 2*sd(mu_star[,3]),  2*sd(mu_star[,4]), cor(mu_star[,3], mu_star[,4])*180, deg = TRUE, border = src.col[2], lwd = 1.5)
		}
	}

}







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# which.eng.setting = 'Setting 1'; which.time.setting = 'Same rates'; sim = 1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
plot_location_scatter = function(which.eng.setting, which.time.setting, sim, method = "BASCS", main = NULL){
	if(is.null(main)) main = paste('Energy:', which.eng.setting, '    and     Time:', which.time.setting)

	f = subset.file(which.eng.setting, which.time.setting)[sim]
	plot(NA,main=main, pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = '', ylab = '')

	if(method == "BASCS"){
		load(gsub('\\_eBASCS', '', f))
		mu_star = t(sapply(main_run[[1]], function(x){x[2:5]}))[sample(1:10000, 2000),]
	}
	if(method == "eBASCS"){
		load(f)
		mu_star = t(sapply(main_run_time[[1]], function(x){x$par[2:5]}))[sample(1:10000, 2000),]
	}

	points(mu_star[,1], mu_star[,2], pch = 19,cex=0.2, col = src.col.alpha[1])
	points(mu_star[,3], mu_star[,4], pch = 19,cex=0.2, col = src.col.alpha[2])

	points(params$positions, pch = 4, cex = 4, lwd = 2)

}





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Draw a contour line for HPD 2D posterior Interval
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# which.eng.setting = 'Setting 4'; which.time.setting = 'Same rates'; method = 'BASCS'

library('KernSmooth')

plot_location_contours = function(which.eng.setting, which.time.setting, method = 'BASCS'){

	fs = subset.file(which.eng.setting, which.time.setting)

	plot(NA,main="", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = '', ylab = '')

	if(method == "BASCS"){
		for(f in fs){
			load(gsub('\\_eBASCS', '', f))
			mu_star_BASCS = t(sapply(main_run[[1]], function(x){x[2:5]}))
			h = apply(mu_star_BASCS, 2, IQR)/10
			# z = kde2d(mu_star_BASCS[,1], mu_star_BASCS[,2], n = 100)
			# contour(z, col = src.col[1], lwd = 1.5, add = TRUE, nlevels = 2)
			# z = kde2d(mu_star_BASCS[,3], mu_star_BASCS[,4], n = 100)
			# contour(z, col = src.col[2], lwd = 1.5, add = TRUE, nlevels = 2)
			z = bkde2D(mu_star_BASCS[,1:2], h[1:2])
			contour(z$x1, z$x2, z$fhat, col = src.col[1], lwd = 1.5, nlevels = 2, labels = '', add = TRUE)
			z = bkde2D(mu_star_BASCS[,3:4], h[3:4])
			contour(z$x1, z$x2, z$fhat, col = src.col[2], lwd = 1.5, nlevels = 2, labels = '', add = TRUE)
		}
	}

	if(method == "eBASCS"){
		for(f in fs){
			load(f)
			mu_star = t(sapply(main_run_time[[1]], function(x){x$par[2:5]}))
			h = apply(mu_star, 2, IQR)/10
			# z = kde2d(mu_star[,1], mu_star[,2], n = 200)
			# contour(z, col = src.col[1], lwd = 1.5, add = TRUE, nlevels = 2)
			# z = kde2d(mu_star[,3], mu_star[,4], n = 200)
			# contour(z, col = src.col[2], lwd = 1.5, add = TRUE, nlevels = 2)
			z = bkde2D(mu_star_BASCS[,1:2], h[1:2])
			contour(z$x1, z$x2, z$fhat, col = src.col[1], lwd = 1.5, nlevels = 2, labels = '', add = TRUE)
			z = bkde2D(mu_star_BASCS[,3:4], h[3:4])
			contour(z$x1, z$x2, z$fhat, col = src.col[2], lwd = 1.5, nlevels = 2, labels = '', add = TRUE)
		}
	}

}




