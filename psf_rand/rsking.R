rsking <- function(draws,r0,slope){
  # NOTE: NO KING.NORM - r0 and slope refer to ENVELOPE ONLY
  
  # Envelope function
  c2d <- function(x){
    dcauchy(x[1],scale=r0)*dcauchy(x[2],scale=r0)
  }
  envelope <- max(king(c(0,0),c(0,0))/c2d(c(0,0)),5)
  count <- 1
  # Rejection sampling
  samples <- matrix(NA,draws,2)
  while (count <= draws){
    proposal <- c(rcauchy(1,scale=r0),rcauchy(1,scale=r0))
    u <- runif(1,0,1)
    if (u < king(proposal,c(0,0))/(envelope*c2d(proposal))){
      samples[count,] <- proposal
      count <- count + 1
    }
  }
  return(samples)
}
