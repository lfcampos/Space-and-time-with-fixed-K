radius <- function(x,y,off.angle,ellip){
  value <- (x*cos(off.angle)+y*sin(off.angle))^2+(y*cos(off.angle)-x*sin(off.angle))^2/(1-ellip)^2
  return(sqrt(value))
}