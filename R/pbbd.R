library(ibd)
psfreq=function(design)
{
  v=max(design)
  k=ncol(design)
  b=nrow(design)
  trtpos.freq=matrix(0,v,k)
  for(pos in 1:k)
  {
    for (blk in 1:b)
    {
      trtpos.freq[design[blk,pos],pos]=trtpos.freq[design[blk,pos],pos]+1
    }
  }
  return(trtpos.freq)
}

cycle = function(x)
{
  x = c(x[2:length(x)],x[1])
  return(x)
}

dpf = function(v,b,r,k)
{
  m = floor(b/v)
  m
  k1 = r-k*m
  k2 = k - k1
  temp = c(rep(m + 1, k1),rep(m, k2))
  dpm = matrix(0,v,k)
  dpm[1,] = temp
  for(i in 2:v) dpm[i,] = cycle(dpm[i-1,])
  return(dpm)
}

allocate = function(j,v,b,k,mvec,x1,x2)
{
  Bs = NULL
  for(i in 1:v)
  {
    m = mvec[i]
    Bi = which(apply(x1[,1:k], 1, function(x) any(x == i)))
    Bia = setdiff(Bi, Bs)
    if(length(Bia) < m) return(0)
    if(length(Bia) == 1) Bs.t = Bia else Bs.t = sample(Bia, size = m)
    x2[Bs.t, j] = i 
    x1[Bs.t,][which(x1[Bs.t,] == i)] = 0
    Bs = union(Bs,Bs.t)
  }
  return(out = list(x1 = x1, x2 = x2))
}  

balancify = function(d1)
{
  v = max(d1)
  b = nrow(d1)
  k = ncol(d1)
  r = b*k/v
  if(r - floor(r) != 0) stop("Design should be equireplicate")
  M = dpf(v,b,r,k)
  d2 = matrix(0,b,k)
  j = 0
  trial = 0
  while(j < k & trial <= 10000)
  {
    trial = trial + 1
    j = j + 1
    mvec = M[,j]
    out = allocate(j,v,b,k,mvec, x1 = d1, x2 = d2)
    if(is.list(out)) 
    {
      d1 = out$x1
      d2 = out$x2
    } else j = j - 1
  }
  if(j == k & trial <= 10000) 
  {
      P = psfreq(d2) 
      result = list(design = d2, P = P)
  }  else result = "Try again"
  return(result)
}

pbbd = function(v,b,k)
{
  if(v > 30 & k > 3) stop("Better to use balancify() function with an input design")
  r = b*k/v
  if(r - floor(r) != 0) stop("Equireplicate design not possible")
  d = ibd(v,b,k)
  if(all(diag(d$conc.mat) != r)) stop("Try again or use balancify() function with an equireplicate design")
  d1 = d$design
  d2 = balancify(d1)
  if(is.list(d2))
  {  
    parameters = c(v = v, b = b, r = r, k = k)
    efficiencies = c(Aeff = d$A.Efficiency, Deff = d$D.Efficiency)
    out = list(parameters = parameters, efficiencies = efficiencies, design = d2$design, P = d2$P)
  } else out = d2  
  return(out)
}