generate1=function(n, d)
{
  x = rnorm(n*d)
  x = 2.5*(pnorm(x)-0.5)
  x = matrix(x, nrow=n, ncol=d, byrow=T)
  m = sin(2*pi*x)
  e = rnorm(n)
#   sigma = sapply(1:n, function(j) (sqrt(d)/2) * ((100 - exp(sum(abs(x[j,]))/d))/(100 + exp(sum(abs(x[j,]))/d))) )
  y = sapply(1:n, function(j) sum(m[j,]) + e[j]) #(sigma[j]*e[j]))
  return(list(y=y, m=m, x=x))
}