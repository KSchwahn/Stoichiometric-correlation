divisors <- function(x){
  #This function is automaticcaly called within the function "ks_stochiometric_correlation".
  #It will give the number of "Blocks", in which the correlation matrix can be divided to speed up the calculation and make parallelization more meaningful.
  a=length(combn(x,2,simplify=F))
  b=length(c(combn(c(1:4),2,simplify = F),combn(c(4:1),2,simplify = F)))
  z <- a*b
  v <- z+x
  y <- seq_len(v)
  #  Modulo division. If remainder is 0 that number is a divisor of x so return it
  out<-y[ v%%y == 0 ]
  return(out[which(out<=x)])
}
