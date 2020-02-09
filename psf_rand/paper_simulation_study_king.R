# Simulate photons from three sources with energy. The energy is simulated in three different ways:
# 1) All different Gammas
# 2) Uniform(0,max_back_energy) background and (potentially) different Gammas for each source
# 3) Uniform(0,max_back_energy) background and the same Gamma for all sources
jid=commandArgs(trailingOnly=T)[1] #this means "take a number from the square brackets from the .bsub file
jid=as.numeric(jid)
 #jid - floor((jid-1)/10)*10
sep.vec <- c(0.5,1,1.5,2)
separation <- sep.vec[4] #(floor((jid-1)/10) + 1)/2


# Fresh start
#rm(list=ls())

# Packages
library("lattice")
library("MASS")
library("mvtnorm")
library("coda")
library("cubature")

setwd("/Users/luis/Dropbox/Disambiguating\ sources\ using\ spatial-time-spectral\ information/Code/Space\ and\ time\ with\ fixed\ K/psf_rand")
#Ssetwd("C:/Users/Dave/Desktop/Overlapping Sources/Paper Simulation/Code")
source("rsking.R")
source("radius.R")
source("king_profile_density.R")
source("king.constant.R")

r0 <- 0.6 # 5.128343
slope <- 1.5 # 1.58749
ellip <- 0.00573533
off.angle <- 0
#increment <- 0.001

for (rel_back in c(0.001,0.01,0.1,1)){
  for (rel_intensity in c(1,2,5,10,50)){
      for (datanum in 1:10){ 

# rel_back = c(0.001,0.01,0.1,1)[1]
# rel_intensity = c(1,2,5,10,50)[1]
# datanum = (1:10)[1]
      # Specify image area
      xup <- 5
      xlow <- -5
      yup <- 5
      ylow <- -5
      max_back_energy <- 10000 # for uniform background
      
      
      # Specify number of sources and either enter positions or randomly generate them
      num <- 2 
      positions <- matrix(NA,num,2)
      positions[1,] <-  c(-separation/2,0)
      positions[2,] <-   c(separation/2,0) 
      
      
      # Specify intensities or randomly generate them
      yl <- (yup - ylow)/2
      xl <- (xup - xlow)/2
      lambda_sim <- matrix(NA,num+1,1)
      lambda_sim[1] <- 1000
      lambda_sim[2] <- lambda_sim[1]/rel_intensity
      lambda_sim[3] <- rel_back*lambda_sim[2]*0.513/10.7   # Cauchy doesn't have variance. Just match mean background to 
      # average number of photons in the region defined by denisty greater than 0.05 (the Cauchy also has very long tails
      # so we don't want to use quantiles).
      
      
      # Mapping from x-axis order to brightness order
      position_bright <- matrix(NA,2,num)
      position_bright[1,] <- c(order(positions[,1]))
      position_bright[2,] <- c(order(lambda_sim[1:num]))
      position_bright <- position_bright[,order(position_bright[1,])]
      rownames(position_bright) <- c("x-axis","dim to bright")
      
      
      # Specify Gamma energy parameters
      energy_paras_v1 <- matrix(NA,2,num+1)
      energy_paras_v1[1,] <- c(3,4,5)
      energy_paras_v1[2,] <- c(0.005,0.004,0.001)
      energy_paras_v2 <- matrix(NA,2,num)
      energy_paras_v2[1,] <- c(3,4) 
      energy_paras_v2[2,] <- c(0.005,0.004) 
      energy_paras_v3 <- c(3,0.005)
      
      #curve(dgamma(x,energy_paras_v1[1,num+1],energy_paras_v1[2,num+1]),xlim=c(0,max_back_energy),xlab="Energy",ylab="Density")
      #x <- seq(0,max_back_energy,1)
      #for (i in 1:num){
      #  points(x,dgamma(x,energy_paras_v1[1,i],energy_paras_v1[2,i]),type="l")
      #}
      
      #curve(dgamma(x,energy_paras_v2[1,1],energy_paras_v2[2,1]),xlim=c(0,max_back_energy),xlab="Energy",ylab="Density")
      #x <- seq(0,max_back_energy,1)
      #points(x,dgamma(x,energy_paras_v2[1,2],energy_paras_v2[2,2]),type="l")
      #y <- rep(1/max_back_energy,length(x))
      #points(x,y,type="l")
      
      #curve(dgamma(x,energy_paras_v3[1],energy_paras_v3[2]),xlim=c(0,max_back_energy),xlab="Energy",ylab="Density")
      #x <- seq(0,max_back_energy,1)
      #y <- rep(1/max_back_energy,length(x))
      #points(x,y,type="l")
      
      
      # Generate from sources
      N_sim <- matrix(NA,num+1,1)
      sall_sim_v1 <- vector("list",num)
      sall_sim_v2 <- vector("list",num)
      sall_sim_v3 <- vector("list",num)
      for (i in 1:num){
        N_sim[i] <- rpois(1,lambda_sim[i])
        king.norm <- 1/king.constant(100)$integral
        idraws <- rsking(N_sim[i],r0,slope)+cbind(rep(positions[i,1],N_sim[i]),rep(positions[i,2],N_sim[i])) 
        sall_sim_v1[[i]] <- cbind(idraws,rgamma(N_sim[i],energy_paras_v1[1,i],energy_paras_v1[2,i]))
        sall_sim_v2[[i]] <- sall_sim_v1[[i]]
        sall_sim_v2[[i]][,3] <- rgamma(N_sim[i],energy_paras_v2[1,i],energy_paras_v2[2,i])
        sall_sim_v3[[i]] <- sall_sim_v1[[i]]
        sall_sim_v3[[i]][,3] <- rgamma(N_sim[i],energy_paras_v3[1],energy_paras_v3[2])
      }
      
      
      # Delete observations outside image
      obs_sim_v1 <- vector("list",num)
      obs_sim_v2 <- vector("list",num)
      obs_sim_v3 <- vector("list",num)
      for (i in 1:num){
        obs_sim_v1[[i]] <- sall_sim_v1[[i]][sall_sim_v1[[i]][,1] >= xlow & sall_sim_v1[[i]][,1] <= xup & sall_sim_v1[[i]][,2] >= ylow & sall_sim_v1[[i]][,2] <= yup,]
        obs_sim_v2[[i]] <- sall_sim_v2[[i]][sall_sim_v2[[i]][,1] >= xlow & sall_sim_v2[[i]][,1] <= xup & sall_sim_v2[[i]][,2] >= ylow & sall_sim_v2[[i]][,2] <= yup,]
        obs_sim_v3[[i]] <- sall_sim_v3[[i]][sall_sim_v3[[i]][,1] >= xlow & sall_sim_v3[[i]][,1] <= xup & sall_sim_v3[[i]][,2] >= ylow & sall_sim_v3[[i]][,2] <= yup,]
      }
      obs_nums <- matrix(NA,num+1,1)
      for (i in 1:num){
        obs_nums[i] <- length(obs_sim_v1[[i]][,1])
      }
      
      #for (i in 1:num){
      #  print(length(sall_sim_v1[[i]][,1]))
      #  print(length(obs_sim_v1[[i]][,1]))
      #  print(length(sall_sim_v2[[i]][,1]))
      #  print(length(obs_sim_v2[[i]][,1]))
      #  print(length(sall_sim_v3[[i]][,1]))
      #  print(length(obs_sim_v3[[i]][,1]))
      #}
      
      # Generate background
      N_sim[num+1] <- rpois(1,4*xl*yl*lambda_sim[num+1])
      xb <- runif(N_sim[num+1],xlow,xup)
      yb <- runif(N_sim[num+1],ylow,yup)
      eb_v1 <- rgamma(N_sim[num+1],energy_paras_v1[1,num+1],energy_paras_v1[2,num+1])
      eb_v2 <- runif(N_sim[num+1],0,max_back_energy)
      eb_v3 <- eb_v2
      obs_nums[num+1] <- N_sim[num+1]
      
      # Collect energies
      energy_v1 <- eb_v1
      for (i in 1:num){
        energy_v1 <- c(energy_v1,obs_sim_v1[[i]][,3])
      }
      energy_v2 <- eb_v2
      for (i in 1:num){
        energy_v2 <- c(energy_v2,obs_sim_v2[[i]][,3])
      }
      energy_v3 <- eb_v3
      for (i in 1:num){
        energy_v3 <- c(energy_v3,obs_sim_v3[[i]][,3])
      }
      
      
      #Plot
      x <- xb
      y <- yb
      for (i in 1:num){
        x <- c(x,obs_sim_v1[[i]][,1])
      }
      for (i in 1:num){
        y <- c(y,obs_sim_v1[[i]][,2])
      }
      pdf(paste("/n/home06/david_jones/Results/Latest/plot_r0",r0,"_slope",slope,"_relb",rel_back,"_relint",rel_intensity,"_sep",separation,"_run",datanum,".pdf",sep=""))
      plot(x,y,pch=19,cex=0.1,main="Two Overlapping Sources",xlim=c(xlow,xup),ylim=c(ylow,yup))
      for (i in 1:num){
        points(positions[i,1],positions[i,2],pch=19,cex=1.5,col=i+1)
      }
      dev.off()
      
      
      # Data
      data <- cbind(x,y)
      
      
      # Save
      save(list=c("data","lambda_sim","N_sim","obs_nums","obs_sim_v1","obs_sim_v2","obs_sim_v3","positions","position_bright","sall_sim_v1","sall_sim_v2","sall_sim_v3","x","xb","xl","xlow","xup","y","yb","yl","ylow","yup","max_back_energy","energy_v1","energy_v2","energy_v3","energy_paras_v1","energy_paras_v2","energy_paras_v3"),file=paste("/n/home06/david_jones/Results/Latest/paper_simulation_r0",r0,"_slope",slope,"_relb",rel_back,"_relint",rel_intensity,"_sep",separation,"_run",datanum,".RData",sep=""))
    }}}