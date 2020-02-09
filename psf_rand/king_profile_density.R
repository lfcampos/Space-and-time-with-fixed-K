# King profile PSF

king <- function(pos,loc){
  m <- length(pos)/2 
  if(m > 1){
  x <- pos[,1]-loc[1]
  y <- pos[,2]-loc[2]
  rad <- radius(x,y,off.angle,ellip)
  value <- king.norm/(1+(rad/r0)^2)^slope
  }
  if(m == 1){
    x <- pos[1]-loc[1]
    y <- pos[2]-loc[2]
    rad <- radius(x,y,off.angle,ellip)
    value <- king.norm/(1+(rad/r0)^2)^slope
  }
  return(value)
}

