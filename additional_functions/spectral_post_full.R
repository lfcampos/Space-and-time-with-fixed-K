spectral_post <- function(m1,a1,m2,a2,ew,dat){
  # Proposals are all symmetric and so can be ignored
  priorm <- -log(emean.range)
  priora <- log(dgamma(a1,ashape,arate))
  value <- sum(log(dgamma(dat,a1,a1/m1)))+priorm+priora
  return(value)
}