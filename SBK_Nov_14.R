#library(corpcor)
library(locpol)
require(corpcor)
source("bspline.constant.R")
source("generate1.R")

z= c(2/3, 1/16)
Mod(polyroot(c(1, -z)))

#number of functions
d=2
id = 1  # which function to estimate
NumNode = c(50, 100, 150, 200)
eff = matrix(nrow = 100, ncol = 4)

generate1=function(n, d)
{
  e = rnorm(n, 0, 1)
  x = e
  m = matrix(nrow=n, ncol = d)
  for(t in (d+1):n)
  {
    m[t,1] = -1/16 * cos(2*pi*x[t-1])
    m[t,2] =  2/(x[t-2]^2+3)#sin(2*pi*x[t-2])#1/(x[t-2]^2+1) - 6/5 #sin(2*pi*x[t-1]) #2/(x[t-1]^2+3)
    x[t] = sum(m[t,])+ e[t]
  }
  
  y = x[(d+1):n]
  x.mat = x[d:(n-1)]
  for(i in 2:d)
  {
    x.mat = cbind(x.mat, x[(d+1-i):(n-i)])
  }
  
  m = m[(d+1):n,]
  
  return(list(y=y, x=x.mat, m=m))
}

for(kk in 1:4){
  
  for(k in 1:100){
    n = NumNode[kk]
    data = generate1(n=n,d=d)
    y = data$y
    x = data$x
    m = data$m
    
    n = length(y)
    
    N = min(floor(n^0.4 * log(n))+1, floor((n/2-1)/d)) #the number of knots
    knots = matrix(nrow=d, ncol = N+2)
    B = matrix(nrow=n, ncol = (N+1)*d)
    for(i in 1:d)
    {
      knots[i,] = seq(min(x[,i]),max(x[,i]), length=N+2)
      B[, ((N+1)*(i-1) + 1):(i*(N+1))] = bspline(x[,i],knots[i,-c(1,N+2)], knots[i,c(1,N+2)])
    }
    lam = pseudoinverse(t(B)%*%B)%*%t(B)%*%y #least square error
    
    # the pre estimated values. m2, m3...
    m.hat = matrix(nrow=n, ncol = d)
    for(i in 1:d)
    {
      m.hat[,i] = as.matrix(B[, ((N+1)*(i-1) + 1):(i*(N+1))])%*%lam[((N+1)*(i-1)+1):(i*(N+1))]
    }
    
    #the backfit stage
    y.p = matrix(nrow=n, ncol = d)
    for(alpha in 1:d)
    {
      
      y.p[,alpha] = y - apply(as.matrix( m.hat[,-alpha] ),1 , sum)
    }
    
    #the oracle
    y.o = matrix(nrow=n, ncol = d)
    
    for(alpha in 1:d)
    {
      y.o[,alpha] = y - apply(as.matrix( m[,-alpha] ),1 , sum)
    }
    mr = m[, id]
    
    
    #The third part, kernel
    bw = .95
    ?thumbBw()
    X = x[,id]
    y = y.p[,id]
    m = y.o[,id]
    x = x[,id]
    
    bw = thumbBw(x,y,1, QuartK)
    bw.o = thumbBw(x,m, 1,QuartK)
    
    m.ll = numeric(n)
    for(j in 1:n)
    {
      xt = X-x[j]
      kalpha = 15/16*((1-(xt/bw)^2+abs(1-(xt/bw)^2))/2)^2/bw
      W=diag(as.numeric(kalpha))
      
      beta = pseudoinverse(t(xt)%*%W%*%xt)%*%t(xt)%*%W%*%y
      
      m.ll[j] = beta[1]
    }
    
    m.o = numeric(n)
    for(j in 1:n)
    {
      xt = X-x[j]
      kalpha = 15/16*((1-(xt/bw.o)^2+abs(1-(xt/bw.o)^2))/2)^2/bw.o
      W=diag(as.numeric(kalpha))
      
      beta = pseudoinverse(t(xt)%*%W%*%xt)%*%t(xt)%*%W%*%m
      m.o[j] = beta[1]
    }
    
    eff[k, kk] = sum((m.o - mr)^2) / sum((m.ll - mr)^2)
    cat("                \r")
    cat(k,"\r")
  }
  cat("                         \r")
  cat(n, "\n")
}

# Kernel Density Plot
par(mfrow = c(2,2))
par(mar=c(2,2,2,2))
globaly = c(0, 2)
re <- density(as.numeric(eff[,1])) # returns the density data 
plot(re, main = "n = 100", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)
re <- density(as.numeric(eff[,2])) # returns the density data 
plot(re, main = "n = 150", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)

re <- density(as.numeric(eff[,3])) # returns the density data 
plot(re, main = "n = 200", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)
re <- density(as.numeric(eff[,4])) # returns the density data 
plot(re, main = "n = 400", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)



