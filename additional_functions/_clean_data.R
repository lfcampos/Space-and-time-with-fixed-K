# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Clean the data 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

spatial <- dat[,1:2]
arrival_time <- dat[,4]
energy <- dat[,3]
max_back_energy <- max(energy)

spatial <- as.matrix(spatial)
colnames(spatial) <- NULL
rownames(spatial) <- NULL

# Basic data information
obs_num <- length(spatial[,1])
xlow <- min(spatial[,1])
xup <- max(spatial[,1])
ylow <- min(spatial[,2])
yup <- max(spatial[,2])
yl <- (yup - ylow)/2
xl <- (xup - xlow)/2
img_area <- 4*xl*yl  # Image area - rectangle assumed
