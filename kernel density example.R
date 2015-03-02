library(KernSmooth)
attach(faithful)
x = waiting
y = eruptions

plot(x, y)
hist(x)

h =5

# rule of thumb bandwidth
h = 1.06*sd(x)*length(x)^(-1/5)

########a close up of the densities

sum.dens = matrix(nrow = length(x), ncol = length(x))
plot(x,rep(0,length(x)), ylim = c(-.01, 1), ylab = "density")




# for the first point plot a density over the point
ind = 1
x1 = x[ind]
x.t = seq(min(x), max(x), length = length(x)*10)
lines(x.t, dnorm(x.t, x1,h)/max(dnorm(x.t, x1,h)), col = "red")
sum.dens[ind,] = dnorm(x, x1,h)


# Now do for the second point
ind = 2
x1 = x[ind]
lines(x.t, dnorm(x.t, x1,h)/max(dnorm(x.t, x1,h)), col = "red")
sum.dens[ind,] = dnorm(x, x1,h)

# keep doing this for all X's
for(ind in 3:length(x))
{
  x1 = x[ind]
  lines(x.t, dnorm(x.t, x1,h)/max(dnorm(x.t, x1,h)), col = "red")
  sum.dens[ind,] = dnorm(x, x1,h)
}




######Now we need to sum the densities to find the the overall density
#we need to plot with larger ylims to see the sum

plot(x,rep(0,length(x)), ylim = c(-.01, max(rowSums(sum.dens))), ylab = "density")


# plot the individual densities
for(ind in 1:length(x))
{
  x1 = x[ind]
  lines(x.t, dnorm(x.t, x1,h), col = "red")
}

#plot the sum
lines(sort(x),rowSums(sum.dens)[order(x)])


