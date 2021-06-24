
# take 5 digit of decimal value and cut the numbers
floorDec <- function(valParm ,x){
  y <- function(x, level=1) round(x - 5*10^(-level-1), level)
  res <-y(as.numeric(valParm),x)
  return(res)
}

