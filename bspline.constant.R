bspline=function(x, inner.knots, Boundary.knots) 
{
  Boundary.knots=sort(Boundary.knots);
  x = as.matrix(x)
  knots=c(Boundary.knots[1], sort(inner.knots),
          Boundary.knots[2]);
  np=length(inner.knots)+1
  s=matrix(nrow = nrow(x), ncol = np)
  for(j in 1:nrow(x))
  {
    if(x[j]==Boundary.knots[2]) 
    {
      s[j,np]=1
      s[j,-np] = 0
    } else {
      for( i in 1: np)
        s[j,i]=basis(x[j], i, knots)
    }
  }
  return(s)
}

basis=function(x, i, knots)
{ 
  
  if((x<knots[i+1])&(x>=knots[i])) y=1 
  else 
    y=0
  return(y)
}