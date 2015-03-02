#library(corpcor)
require(corpcor)
source("bspline.constant.R")
source("generate1.R")

d=10
id = 2  # which function to estimate
NumNode = c(50, 50, 50, 50)
eff = matrix(nrow = 50, ncol = 4)

for(kk in 1:4){
  n = NumNode[kk]
for(k in 1:50){

# id is which x that we want to estimate eventually
data = generate1(n=n,d=d)
y = data$y
x = data$x
m = data$m

#the pre-estimate stage (splines)
#N=20   #!!!! use the equation in paper 2486
N = min(floor(n^0.4 * log(n))+1, floor((n/2-1)/d))
knots = matrix(nrow=d, ncol = N+2)
B = matrix(nrow=n, ncol = (N+1)*d)
for(i in 1:d)
{
  knots[i,] = seq(min(x[,i]),max(x[,i]), length=N+2)
  B[, ((N+1)*(i-1) + 1):(i*(N+1))] = bspline(x[,i],knots[i,-c(1,N+2)], knots[i,c(1,N+2)])
}

lam = pseudoinverse(t(B)%*%B)%*%t(B)%*%y #least square error

m.hat = matrix(nrow=n, ncol = d)
for(i in 1:d)
{
  m.hat[,i] = B[, ((N+1)*(i-1) + 1):(i*(N+1))]%*%lam[((N+1)*(i-1)+1):(i*(N+1))]
}
# the pre estimated values. m2, m3

#the backfit stage
y.p = matrix(nrow=n, ncol = d)
for(alpha in 1:d)
{
  y.p[,alpha] = y - apply(as.matrix(m.hat[,-alpha]),1 , sum)
}


#the oracle
y.o = matrix(nrow=n, ncol = d)
for(alpha in 1:d)
{
  y.o[,alpha] = y - apply(as.matrix(m[,-alpha]),1 , sum)
}

#plot(x[,id], y, col= "grey")
#plot true funciton
#lines(sort(x[,id]), m[order(x[,id]),id])
mr = m[order(x[,id]), id]
#plot preestimate
#lines(sort(x[,id]), m.hat[order(x[,id]),id], col="red")
#plot pseudo-response (backfit)
#points(sort(x[,id]), y.p[order(x[,id]),id], col="blue")
#plot oracle
#points(sort(x[,id]), y.o[order(x[,id]),id], col="green")

########### connect two parts

#SumCol = sapply(1:n, function(i) sum(m.hat[i, ]))
#SumCol = SumCol - m.hat[,id]
#y.p = y.p - SumCol
###########

bw = .95
X = sort(x[,id])
y = y.p[order(x[,id]),id]
m = y.o[order(x[,id]),id]
x = sort(x[,id])

m.ll = numeric(n)
for(j in 1:n)
{
  xt = X-x[j]
  kalpha = 15/16*((1-(xt/bw)^2+abs(1-(xt/bw)^2))/2)^2/bw
  W=diag(as.numeric(kalpha))
  
  beta = pseudoinverse(t(xt)%*%W%*%xt)%*%t(xt)%*%W%*%y
  m.ll[j] = beta[1]
}
#m.ll
m.o = numeric(n)
for(j in 1:n)
{
  xt = X-x[j]
  kalpha = 15/16*((1-(xt/bw)^2+abs(1-(xt/bw)^2))/2)^2/bw
  W=diag(as.numeric(kalpha))
  
  beta = pseudoinverse(t(xt)%*%W%*%xt)%*%t(xt)%*%W%*%m
  m.o[j] = beta[1]
}


eff[k, kk] = sum((m.o - mr)^2) / sum((m.ll - mr)^2)
cat("                        \r")
cat(k,"\r")
}
cat("                         \r")
cat(n, "\n")
}
# Kernel Density Plot
par(mfrow = c(2,2))
par(mar=c(2,2,2,2))
globaly = c(0, 3)
re <- density(as.numeric(eff[,1])) # returns the density data 
plot(re, main = "n = 100", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)
re <- density(as.numeric(eff[,2])) # returns the density data 
plot(re, main = "n = 200", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)

re <- density(as.numeric(eff[,3])) # returns the density data 
plot(re, main = "n = 500", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)
re <- density(as.numeric(eff[,4])) # returns the density data 
plot(re, main = "n = 1000", xlim = globaly) # plots the results
abline(v = 1, col=3.5, lty = 2)


