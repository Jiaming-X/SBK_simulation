#Question 1
#version 1
density = matrix(nrow = 3, ncol = 2)
density[,1] = 0
density[,2] = 1
p = c(0.3, 0.2, 0.5)

MixSample<-function(n = n, p = p, density = density){
  k = length(p)
  y = sample(1:k, n, replace = TRUE, prob = p)
  y2 = matrix(nrow = n, ncol = 1)
  for(i in 1:n){
    y2[i,1] = rnorm(n=1, mean=density[y[i],1], sd=density[y[i],2])
  }
  return(y2)
}
temp = matrix(nrow = 100, ncol = 3)
for(i in 1:100){
temp[i, ] = MixSample(n = 3, p = p, density = density)
}

#version 2
#Three distributions; 
density = matrix(nrow = 3, ncol = 2)
density[,1] = 0
density[,2] = 1
p = c(0.3, 0.2, 0.5)

MixSample<-function(n = n, p = p, density = density){
  k = length(p)
  y = sample(1:k, n, replace = TRUE, prob = p)
  y2 = matrix(nrow = n, ncol = 1)
  for(i in 1:n){
    y2[i,1] = rnorm(n=1, mean=density[y[i],1], sd=density[y[i],2])
  }
  return(y2)
}
temp = matrix(nrow = 100, ncol = 3)
for(i in 1:100){
  temp[i, ] = MixSample(n = 3, p = p, density = density)
}

