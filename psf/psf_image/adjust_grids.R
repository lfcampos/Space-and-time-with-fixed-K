adjust_grids <- function(){
  
  lenx <- length(xpsf_grid)+3
  xpsf_grid2 <- matrix(NA,lenx,1)
  xpsf_grid2[1] <- -10^10
  xpsf_grid2[lenx] <- 10^10
  xpsf_grid2[2] <- xpsf_grid[1] - 10^(-6)
  xpsf_grid2[lenx-1] <- xpsf_grid[lenx-3] + 10^(-6)
  hx <- psf_image[[1]][2]-psf_image[[1]][1]
  xpsf_grid2[3:(lenx-2)] <- xpsf_grid[-(lenx-3)]+hx/2
  
  leny <- length(ypsf_grid)+3
  ypsf_grid2 <- matrix(NA,leny,1)
  ypsf_grid2[1] <- -10^10
  ypsf_grid2[leny] <- 10^10
  ypsf_grid2[2] <- ypsf_grid[1] - 10^(-6)
  ypsf_grid2[leny-1] <- ypsf_grid[leny-3] + 10^(-6)
  hy <- psf_image[[2]][2]-psf_image[[2]][1]
  ypsf_grid2[3:(leny-2)] <- ypsf_grid[-(leny-3)]+hy/2
  
  psf_image_values2 <- matrix(10^(-50),lenx-1,leny-1) # Need to avoid density values of exactly zero
  psf_image_values2[-c(1,lenx-1),-c(1,leny-1)] <- psf_image_values
  # Main PSF
  norm_con <- sum(psf_image_values*hx*hy) 
  # Very low density perimeter (contributes negligible density)
  norm_con <- norm_con + 4*(xup-xlow-hx*(lenx-3)/2)*(yup-ylow-hy*(leny-3)/2)*10^(-50) 
  norm_con <- norm_con + 2*hx*(lenx-3)*(yup-ylow-hy*(leny-3)/2)*10^(-50) 
  norm_con <- norm_con + 2*hy*(leny-3)*(xup-xlow-hx*(lenx-3)/2)*10^(-50)
  psf_image_values2 <- psf_image_values2/norm_con

  value <- list(xpsf_grid2,ypsf_grid2,psf_image_values2)
  return(value)
}