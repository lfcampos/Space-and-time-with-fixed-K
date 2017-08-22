psf2 <- function(pts,center){
  
  m <- length(pts)/2 
  x <- pts[1:m]-center[1]
  y <- pts[(m+1):(2*m)]-center[2]
  grid_indexes <- cbind(findInterval(x,vec=xpsf_grid),findInterval(y,vec=ypsf_grid))
  value <- 0 
  if (m > 1){
    interpx1 <- xpsf_grid[grid_indexes[,1]]
    interpx2 <- xpsf_grid[grid_indexes[,1]+1]
    interpy1 <- ypsf_grid[grid_indexes[,2]]
    interpy2 <- ypsf_grid[grid_indexes[,2]+1]
    c_common <- (interpx2-interpx1)*(interpy2-interpy1)
    c1 <- (interpx2-x)*(interpy2-y)
    c2 <- (x-interpx1)*(interpy2-y)
    c3 <- (interpx2-x)*(y-interpy1)
    c4 <- (x-interpx1)*(y-interpy1)
    value <- (1/c_common)*(c1*psf_image_values[grid_indexes]+c2*psf_image_values[cbind(grid_indexes[,1]+1,grid_indexes[,2])]+c3*psf_image_values[cbind(grid_indexes[,1],grid_indexes[,2]+1)]+c4*psf_image_values[grid_indexes+1])
  } 
  if (m == 1){
    interpx1 <- xpsf_grid[grid_indexes[1]]
    interpx2 <- xpsf_grid[grid_indexes[1]+1]
    interpy1 <- ypsf_grid[grid_indexes[2]]
    interpy2 <- ypsf_grid[grid_indexes[2]+1]
    c_common <- (interpx2-interpx1)*(interpy2-interpy1)
    value <- (1/c_common)*(c1*psf_image_values[grid_indexes]+c2*psf_image_values[cbind(grid_indexes[1]+1,grid_indexes[2])]+c3*psf_image_values[cbind(grid_indexes[1],grid_indexes[2]+1)]+c4*psf_image_values[grid_indexes+1])
  }

  return(value)
}