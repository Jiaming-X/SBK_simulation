require(corpcor)
source("bspline.constant.R")
source("generate1.R")

n=400
d=10
data = generate1(n=n,d=d)

y = data$y
x = data$x
m = data$m

#the pre-estimate stage (splines)
N=20
knots = matrix(nrow=d, ncol = N+2)
B = matrix(nrow=n, ncol = (N+1)*d)
for(i in 1:d)
{
  knots[i,] = seq(min(x[,i]),max(x[,i]), length=N+2)
  B[, ((N+1)*(i-1) + 1):(i*(N+1))] = bspline(x[,i],knots[i,-c(1,N+2)], knots[i,c(1,N+2)])
}

lam = pseudoinverse(t(B)%*%B)%*%t(B)%*%y

m.hat = matrix(nrow=n, ncol = d)
for(i in 1:d)
{
  m.hat[,i] = B[, ((N+1)*(i-1) + 1):(i*(N+1))]%*%lam[((N+1)*(i-1)+1):(i*(N+1))]
}


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

id = 2
plot(x[,id], y, col= "grey")

#plot true funciton
lines(sort(x[,id]), m[order(x[,id]),id])

#plot preestimate
lines(sort(x[,id]), m.hat[order(x[,id]),id], col="red")

#plot pseudo-response (backfit)
points(sort(x[,id]), y.p[order(x[,id]),id], col="blue")

#plot oracle
points(sort(x[,id]), y.o[order(x[,id]),id], col="green")


