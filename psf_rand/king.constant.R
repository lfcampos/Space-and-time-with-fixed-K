king.constant <- function(bound){
  #library("cubature")
  king_fun <- function(pos){
    x <- pos[1]
    y <- pos[2]
    rad <- radius(x,y,off.angle,ellip)
    value <- 1/(1+(rad/r0)^2)^slope
    return(value)
  }
  value <- adaptIntegrate(king_fun,c(-bound,-bound),c(bound,bound))
  return(value)
}

#bound <- 10
#inc <- 0.1
#x2 <- seq(-bound,bound,inc)
#y2 <- seq(-bound,bound,inc)
#z2 <- outer(x2,y2,king2)
#1/sum(z2*inc^2)
