divisors <- function(x){
  #  Vector of numberes to test against
  a=length(combn(x,2,simplify=F))
  b=length(c(combn(c(1:4),2,simplify = F),combn(c(4:1),2,simplify = F)))
  z <- a*b
  v <- z+x
  y <- seq_len(v)
  #  Modulo division. If remainder is 0 that number is a divisor of x so return it
  out<-y[ v%%y == 0 ]
  return(out[which(out<=x)])
}
