bspline=function(x, degree, inner.knots, Boundary.knots) 
{
  Boundary.knots=sort(Boundary.knots);
  x = as.matrix(x)
  knots=c(rep(Boundary.knots[1], (degree+1)), sort(inner.knots),
          rep(Boundary.knots[2], (degree+1)));
  # sort the knots.
  np=degree+length(inner.knots)+1
  s=matrix(nrow = nrow(x), ncol = np)
  for(j in 1:nrow(x))
  {
    if(x[j]==Boundary.knots[2]) 
    {
      s[j,np]=1
      s[j,-np] = 0
    } else {
      for( i in 1: np)
        s[j,i]=basis(x[j], degree, i, knots)
    }
  }
  return(s)
}

basis=function(x, degree, i, knots)
{ 
  if(degree==0)
  { 
    if((x<knots[i+1])&(x>=knots[i])) y=1 
    else 
      y=0
  }else{     # how to calculate the value of y ???
    if((knots[degree+i]-knots[i])==0) {
      temp1=0
    } else {
      temp1=(x-knots[i])/(knots[degree+i]-knots[i])
    };
    if((knots[i+degree+1]-knots[i+1])==0) 
    {
      temp2=0
    } else {
      temp2=(knots[i+degree+1]-x)/(knots[i+degree+1]-knots[i+1])
    }
    y= temp1*basis(x, (degree-1), i, knots) +temp2*basis(x, (degree-1),(i+1), knots)
  }
  return(y)
}


