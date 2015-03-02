library(SemiPar)
data(lidar, package="SemiPar")
x <- lidar$range
y = lidar$logratio
plot(x,y, col = "grey", main = "Lidar Degree 1 Spline")

N =20


knots = seq(min(x),max(x),length = N+2)
# knots = c(390, seq(550, max(x), length = 20))
# N = length(knots) - 2
inner.knots = knots[c(-1,-(N+2))]
Boundary.knots = knots[c(1,(N+2))]

#generate the linear b-spline matrix
B = bspline(x, 3, inner.knots, Boundary.knots)

library(corpcor)
b.hat = pseudoinverse(t(B)%*%B)%*%t(B)%*%y
f.hat = B%*%b.hat

lines(x,f.hat, col = "red")


