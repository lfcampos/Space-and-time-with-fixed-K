psf <- function(pts,center){
  m <- length(pts)/2 
  x <- pts[1:m]-center[1]
  y <- pts[(m+1):(2*m)]-center[2]
  grid_indexes <- cbind(findInterval(x,vec=xpsf_grid),findInterval(y,vec=ypsf_grid))
  value <- psf_image_values[grid_indexes]
  return(value)
}