library(corpcor)

n = 100
d = 1
data = generate1(n, d)

bw = .95
X = data$x
y = data$y
m = data$m

x = data$x


m.ll = numeric(n)
for(j in 1:n)
{
  xt = X-x[j]
  kalpha = 15/16*((1-(xt/bw)^2+abs(1-(xt/bw)^2))/2)^2/bw
  W=diag(as.numeric(kalpha))
  
  beta = pseudoinverse(t(xt)%*%W%*%xt)%*%t(xt)%*%W%*%y
  m.ll[j] = beta[1]
}

