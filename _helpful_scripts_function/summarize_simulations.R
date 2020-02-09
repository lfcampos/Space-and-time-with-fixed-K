

file.names = as.character(read.table('file.names.txt')[,1])
file.names = paste('./results/', file.names, sep = '')


figures.names = gsub('results', 'sim_analysis', file.names)
figures.names = gsub('Rda', 'pdf', figures.names)

which.eBASCS = grep('eBASCS', file.names)

file.eBASCS = file.names[ which.eBASCS]
file.BASCS  = file.names[-which.eBASCS]

figures.eBASCS = figures.names[ which.eBASCS]
figures.BASCS  = figures.names[-which.eBASCS]



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Go through and make individual figures for results (same as before)
# More importantly, keep some output and make some joint plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
library(RColorBrewer)
# 
# bright, dim, background
src.col =  c(brewer.pal(9, 'Set1')[1:2], "#000000")
src.col.alpha = paste(src.col, c('4D', '4D', '1A'), sep = '') # add alpha value 0.3, 
src.lab = c('Bright Source', 'Dim Source', 'Background')


source('./additional_functions/light_curves.R')
source('./summarize_functions.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# label Simulations by Setting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

eng.setting = c('equalE','varyE_bySource','varyE_byTime_Source','varyE_byTime_Source_same_Marginal')

file.params = list()

file.params$name = file.eBASCS
file.params$eng.setting = rep('', length(file.eBASCS))
file.params$eng.setting[grep(eng.setting[1], file.eBASCS)] = 'Setting 1'
file.params$eng.setting[grep(eng.setting[2], file.eBASCS)] = 'Setting 2'
file.params$eng.setting[grep(eng.setting[3], file.eBASCS)] = 'Setting 3'
file.params$eng.setting[grep(eng.setting[4], file.eBASCS)] = 'Setting 4'

file.params$sim_number = as.numeric(gsub('\\./results/analysis\\_eBASCS\\_', '', gsub('\\_param_set\\_|\\.Rda', '', gsub(eng.setting[1], '', gsub(eng.setting[2], '', gsub(eng.setting[3], '', gsub(eng.setting[4], '', file.eBASCS)))))))

file.params$time.setting = rep('', length(file.eBASCS))
file.params$time.setting[sim_number %in% 1:10] =  'Same rates'
file.params$time.setting[sim_number %in% 11:20] = 'Diff rates'




home = '.'

function_dir = paste(home, '/additional_functions/', sep = '')
source(paste(function_dir, 'mcmc.spatial.time.R', sep = ''))
source(paste(function_dir, 'mcmc.spatial.R', sep = ''))
source(paste(function_dir, 'light_curves.R', sep = ''))
source(paste(function_dir, "log_posterior_spatial.R",sep=""))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Load Run Analysis Function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source(paste(function_dir, 'run_analysis_summary.R', sep = ''))


n.sims = length(file.params$name)

sim_params = vector('list', n.sims)
sim_analysis_time = vector('list', n.sims)
sim_analysis = vector('list', n.sims)

for(sim in 1:length(file.params$name)){
	print(paste('Analysis ', sim, ': ', file.params$name[sim], sep = ''))
	sim_params[[sim]] = params
	load(file = file.params$name[sim])
	sim_analysis_time[[sim]] = run_analysis_summary_eBASCS(dat, main_run_time, params = params, inc.time = TRUE, make.fig = FALSE, to.pdf = TRUE)
}

for(sim in 1:length(file.params$name)){
	fn = gsub('_eBASCS', '', file.params$name[sim])
	print(paste('Analysis ', sim, ': ', fn, sep = ''))
	sim_params[[sim]] = params
	load(file = fn)
	sim_analysis[[sim]] = run_analysis_summary(dat, main_run, params = params, inc.time = FALSE, make.fig = FALSE, to.pdf = TRUE)
}


save(list = c("sim_params","sim_analysis_time","sim_analysis"), file = './sim_analysis/simulation_summary.Rda')


load(file = './sim_analysis/simulation_summary.Rda')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Make figures only
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
counter = 1
eng.setting = c('Setting_1_A', 'Setting_1_B', 'Setting_2_A', 'Setting_2_B',
	'Setting_3_A', 'Setting_3_B', 'Setting_4_A', 'Setting_4_B')

for(sim in c(1, 3, 21, 23, 41, 43, 61, 63)){
	print(paste('Analysis ', sim, ': ', file.params$name[sim], sep = ''))
	load(file = file.params$name[sim])

	pdf(file = paste('./sim_analysis/analysis_eBASCS_', eng.setting[counter], '.pdf', sep = ''))
	tmp = run_analysis_summary_eBASCS(dat, main_run_time, params = params, inc.time = TRUE, make.fig = TRUE, to.pdf = TRUE)
	dev.off()

	load(file = gsub('_eBASCS', '', file.params$name[sim]))
	pdf(file = paste('./sim_analysis/analysis_BASCS_', eng.setting[counter], '.pdf', sep = ''))
	tmp = run_analysis_summary(dat, main_run, params = params, inc.time = FALSE, make.fig = TRUE, to.pdf = TRUE)
	dev.off()
	counter = counter + 1
}









# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Confidence band plotting function wrapper
# very ugly, need to make cleaner, add legend
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
require(plotrix)
plot.post_CI = function(x_loc, true_vals, post_int, xlab = NULL, ylab = NULL, main = NULL, ylim = NULL, true.unique = FALSE){
	plotCI(x_loc, post_int[1,], ui=post_int[3,], li=post_int[2,], xaxt="n", xlab = xlab, ylab = ylab, main = main, ylim = ylim)
	if(!is.null(true_vals)){
		if(!true.unique) segments(x_tic-5*0.06, true_vals, x_tic+4*0.06, true_vals, col = 'blue', lwd = 2)
		if(true.unique) points(x_loc, true_vals, pch = 4, col = rgb(0, 0, 1, alpha = 0.6))
	}
	axis(1, at=x_tic,labels=x_labels)
	points(x_loc, post_int[1,], pch = 19, col = rgb(1, 0, 0, alpha = 0.6))

	return(NULL)
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Summarize things accross data replicates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

x_groups = do.call('c', lapply(1:4, rep, 10))
x_loc = x_groups + rep((-5:4)*0.06, 4)
x_tic = 1:4


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Location
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

pdf(file = './sim_analysis/_analysis_simulation_Separation.pdf', width = 10)


for(time_setting in c('Same rates', 'Diff rates')){
	subset = which(file.params$time.setting == time_setting)

	dist_sim = sapply(sim_params[subset], function(x) x$separation)
	dist_results_time = sapply(sim_analysis_time[subset], function(x){ c(x$post_mean_dist, x$post_CI_dist) })
	dist_results = sapply(sim_analysis[subset], function(x){ c(x$post_mean_dist, x$post_CI_dist) })
	true.sep = rep(1, 4)
	x_labels = paste('Setting', 1:4)

	ylim = c(0, 2)
	par(mfrow = c(1, 2), oma=c(0,0,2,0)); 
	plot.post_CI(x_loc, true.sep, dist_results_time, xlab = 'Separation', ylab = 'Posterior: Separation', main = 'eBASCS', ylim = ylim, true.unique = FALSE)
	plot.post_CI(x_loc, true.sep, dist_results, xlab = 'Separation', ylab = 'Posterior: Separation', main = 'BASCS', ylim = ylim, true.unique = FALSE)
	title(paste("Posterior CI: Separation -", time_setting), outer=TRUE)
}


dev.off()


for(time_setting in c('Same rates', 'Diff rates')){
	subset = which(file.params$time.setting == time_setting)

	w_sim = sapply(sim_params[subset], function(x) x$obs_nums/sum(x$obs_nums))

	# Space and Time:
	#   a bit more complex, 3 per simulation
	w_results_mean_time = sapply(sim_analysis_time[subset], function(x){ x$post_mean_w})
	w_resultsCI_time = sapply(sim_analysis_time[subset], function(x){ c(x$post_CI_w[,1], x$post_CI_w[,2], x$post_CI_w[,3])})
	w_results_time = list()
	w_results_time[[1]] = rbind(w_results_mean_time[1,], w_resultsCI_time[1:2,])
	w_results_time[[2]] = rbind(w_results_mean_time[2,], w_resultsCI_time[3:4,])
	w_results_time[[3]] = rbind(w_results_mean_time[3,], w_resultsCI_time[5:6,])

	# Space Only:
	w_results_mean = sapply(sim_analysis[subset], function(x){ x$post_mean_w})
	w_resultsCI = sapply(sim_analysis[subset], function(x){ c(x$post_CI_w[,1], x$post_CI_w[,2], x$post_CI_w[,3])})
	w_results = list()
	w_results[[1]] = rbind(w_results_mean[1,], w_resultsCI[1:2,])
	w_results[[2]] = rbind(w_results_mean[2,], w_resultsCI[3:4,])
	w_results[[3]] = rbind(w_results_mean[3,], w_resultsCI[5:6,])

	pdf(file = paste('./sim_analysis/_analysis_simulation_weight', time_setting, '.pdf', sep = ''), width = 10)

	for(i in 1:3){
		par(mfrow = c(1, 2), oma=c(0,0,2,0)); 
		plot.post_CI(x_loc, w_sim[i,], w_results_time[[i]], xlab = 'Separation', ylab = 'Posterior: Source Weight', main = 'eBASCS', ylim = c(0, 1), true.unique = TRUE)
		plot.post_CI(x_loc, w_sim[i,], w_results[[i]], xlab = 'Separation', ylab = 'Posterior: Source Weight', main = 'BASCS', ylim = c(0, 1), true.unique = TRUE)
		title(paste("Posterior CI: Weight -", src.lab[i], '-', time_setting), outer=TRUE)
	}

	dev.off()

}






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Summarize Location
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

for(k in 1:length(file.eBASCS)){
	load(file.eBASCS[k])

	mu_star = t(sapply(main_run_time[[1]], function(x){x$par[2:5]}))

	load(gsub('\\_eBASCS', '', file.eBASCS[k]))

	mu_star_BASCS = t(sapply(main_run[[1]], function(x){x[2:5]}))

	par(mfrow = c(1, 2))

	plot(NA,main="Posterior: Source centers", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')
	points(mu_star[,1], mu_star[,2], pch = 19,cex=0.2, col = src.col[1])
	points(mu_star[,3], mu_star[,4], pch = 19,cex=0.2, col = src.col[2])


	plot(NA,main="Posterior: Source centers", pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = 'x', ylab = 'y')
	points(mu_star_BASCS[,1], mu_star_BASCS[,2], pch = 19,cex=0.2, col = src.col[1])
	points(mu_star_BASCS[,3], mu_star_BASCS[,4], pch = 19,cex=0.2, col = src.col[2])

	readline()
}








# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Posterior Draws plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

pdf(file = 'Location Posterior Results_BASCS_Same_Rates.pdf', height = 7, width = 14)

	for(k in 1:5){
		# png(file = paste('Location Posterior Results_BASCS_', k, '.png', sep = ''), width = 2*480, height = 1*480)
		par(mfrow = c(2, 4), mar=c(3.1,2.1,1.5,0.5), oma = c(0,0,0,0))

		plot_location_scatter('Setting 1', 'Same rates', k, method = "BASCS", main = 'Setting 1')
		plot_location_scatter('Setting 2', 'Same rates', k, method = "BASCS", main = 'Setting 2')
		plot_location_scatter('Setting 3', 'Same rates', k, method = "BASCS", main = 'Setting 3')
		plot_location_scatter('Setting 4', 'Same rates', k, method = "BASCS", main = 'Setting 4')
		plot_location_scatter('Setting 1', 'Same rates', k, method = "eBASCS", main = 'Setting 1')
		plot_location_scatter('Setting 2', 'Same rates', k, method = "eBASCS", main = 'Setting 2')
		plot_location_scatter('Setting 3', 'Same rates', k, method = "eBASCS", main = 'Setting 3')
		plot_location_scatter('Setting 4', 'Same rates', k, method = "eBASCS", main = 'Setting 4')
		# dev.off()
	}
dev.off()

pdf(file = 'Location Posterior Results_eBASCS_Diff_Rates.pdf', height = 7, width = 14)

	for(k in 1:5){
		# png(file = paste('Location Posterior Results_eBASCS_', k, '.png', sep = ''), width = 2*480, height = 1*480)

		par(mfrow = c(2, 4), mar=c(3.1,2.1,1.5,0.5), oma = c(0,0,0,0))

		plot_location_scatter('Setting 1', 'Diff rates', k, method = "BASCS", main = '')
		plot_location_scatter('Setting 2', 'Diff rates', k, method = "BASCS", main = '')
		plot_location_scatter('Setting 3', 'Diff rates', k, method = "BASCS", main = '')
		plot_location_scatter('Setting 4', 'Diff rates', k, method = "BASCS", main = '')

		plot_location_scatter('Setting 1', 'Diff rates', k, method = "eBASCS", main = '')
		plot_location_scatter('Setting 2', 'Diff rates', k, method = "eBASCS", main = '')
		plot_location_scatter('Setting 3', 'Diff rates', k, method = "eBASCS", main = '')
		plot_location_scatter('Setting 4', 'Diff rates', k, method = "eBASCS", main = '')
		# dev.off()
	}
dev.off()




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Posterior Ellipses view 2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


pdf(file = 'Location Posterior Results_BB.pdf', height = 7, width = 10)
	plot_location_ellipses('Setting 1', 'Same rates')
	plot_location_ellipses('Setting 2', 'Same rates')
	plot_location_ellipses('Setting 3', 'Same rates')
	plot_location_ellipses('Setting 4', 'Same rates')

	plot_location_ellipses('Setting 1', 'Diff rates')
	plot_location_ellipses('Setting 2', 'Diff rates')
	plot_location_ellipses('Setting 3', 'Diff rates')
	plot_location_ellipses('Setting 4', 'Diff rates')
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Posterior Ellipses view 2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

pdf(file = 'Location Posterior Results_view2.pdf', height = 7, width = 14)

par(mfrow = c(2, 4), mar=c(3.1,2.1,1.5,0.5), oma = c(0,0,0,0))

	plot_location_ellipses2('Setting 1', 'Same rates', method = 'BASCS', main = 'Setting 1')
	plot_location_ellipses2('Setting 2', 'Same rates', method = 'BASCS', main = 'Setting 2')
	plot_location_ellipses2('Setting 3', 'Same rates', method = 'BASCS', main = 'Setting 3')
	plot_location_ellipses2('Setting 4', 'Same rates', method = 'BASCS', main = 'Setting 4')

	plot_location_ellipses2('Setting 1', 'Diff rates', method = 'BASCS')
	plot_location_ellipses2('Setting 2', 'Diff rates', method = 'BASCS')
	plot_location_ellipses2('Setting 3', 'Diff rates', method = 'BASCS')
	plot_location_ellipses2('Setting 4', 'Diff rates', method = 'BASCS')


par(mfrow = c(2, 4), mar=c(3.1,2.1,1.5,0.5), oma = c(0,0,0,0))

	plot_location_ellipses2('Setting 1', 'Same rates', method = 'eBASCS', main = 'Setting 1')
	plot_location_ellipses2('Setting 2', 'Same rates', method = 'eBASCS', main = 'Setting 2')
	plot_location_ellipses2('Setting 3', 'Same rates', method = 'eBASCS', main = 'Setting 3')
	plot_location_ellipses2('Setting 4', 'Same rates', method = 'eBASCS', main = 'Setting 4')

	plot_location_ellipses2('Setting 1', 'Diff rates', method = 'eBASCS')
	plot_location_ellipses2('Setting 2', 'Diff rates', method = 'eBASCS')
	plot_location_ellipses2('Setting 3', 'Diff rates', method = 'eBASCS')
	plot_location_ellipses2('Setting 4', 'Diff rates', method = 'eBASCS')


dev.off()








# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Posterior Contours
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

pdf(file = 'Location Posterior Results_view2.pdf', height = 7, width = 14)

par(mfrow = c(2, 4), mar=c(3.1,2.1,0.5,0.5), oma = c(0,0,3.5,0))

	plot_location_contours('Setting 1', 'Same rates', method = 'BASCS')
	plot_location_contours('Setting 2', 'Same rates', method = 'BASCS')
	plot_location_contours('Setting 3', 'Same rates', method = 'BASCS')
	plot_location_contours('Setting 4', 'Same rates', method = 'BASCS')

	plot_location_contours('Setting 1', 'Diff rates', method = 'BASCS')
	plot_location_contours('Setting 2', 'Diff rates', method = 'BASCS')
	plot_location_contours('Setting 3', 'Diff rates', method = 'BASCS')
	plot_location_contours('Setting 4', 'Diff rates', method = 'BASCS')


par(mfrow = c(2, 4), mar=c(3.1,2.1,0.5,0.5), oma = c(0,0,3.5,0))

	plot_location_contours('Setting 1', 'Same rates', method = 'eBASCS')
	plot_location_contours('Setting 2', 'Same rates', method = 'eBASCS')
	plot_location_contours('Setting 3', 'Same rates', method = 'eBASCS')
	plot_location_contours('Setting 4', 'Same rates', method = 'eBASCS')

	plot_location_contours('Setting 1', 'Diff rates', method = 'eBASCS')
	plot_location_contours('Setting 2', 'Diff rates', method = 'eBASCS')
	plot_location_contours('Setting 3', 'Diff rates', method = 'eBASCS')
	plot_location_contours('Setting 4', 'Diff rates', method = 'eBASCS')


dev.off()




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Summarize Allocations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

which.eng.setting = 'Setting 4'; which.time.setting = 'Same rates'; method = 'eBASCS'


	fs = subset.file(which.eng.setting, which.time.setting)

	# plot(NA,main=main, pch=19,cex=0.1, xlim = params$positions[,1]*2, ylim = params$positions[,1]*2, xlab = '', ylab = '')

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


			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
			# Summarize Allocation for BASCS
			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
			confusion_star = lapply(main_run_time[[1]], function(x){x$confusion})
			true_alloc = apply(confusion_star[[3]], 2, sum)
			denom = matrix(rep(true_alloc, 3), byrow = TRUE, ncol = 3)
			correct_alloc = sapply(confusion_star, function(x) diag(x/denom))
			par(mfrow = c(1, 3))
			hist(correct_alloc[1,])
			hist(correct_alloc[2,])
			hist(correct_alloc[3,])

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
			# Summarize posterior estimate of intensity
			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
			pi_star = do.call('rbind', lapply(main_run_time[[1]], function(x){x$par[6:8]}))
			apply(pi_star, 2, mean)			

		}
	}



