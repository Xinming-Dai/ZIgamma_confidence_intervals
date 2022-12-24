#This is the simulation of Fiducial method.
# ci.f function is the confidence intverval function
source('functions/zigamma_data.R')

# generalized Fiducial quantities
GFQzig = function(ni0, ni1, zigdata1, nr){
  # GFQ of a
  y=zigdata1^(1/3)
  u=mean(y)
  s2=var(y)
  Z=rnorm(nr,0,1);k=(ni1-1)/rchisq(nr,df=ni1-1)
  Gmu=u+Z*sqrt(s2/ni1)*sqrt(k)
  Gsigma2=s2*k
  Ga=1/9*((1+.5*Gmu^2/Gsigma2)+((1+.5*Gmu^2/Gsigma2)^2-1)^(1/2))
  
  # GFQ of delta
  n=ni0+ni1
  ZW=rnorm(nr,0,1); ZW2=ZW^2
  Gdelta=(ni0+ZW2/2)/(n+ZW2)-ZW*sqrt(ni0*(1-ni0/n)+ZW2/4)/(n+ZW2)
  return(list(Gmu=Gmu,Gsigma2=Gsigma2,Ga=Ga,Gdelta=Gdelta))
}

# generalized Fiducial quantitie of Ggamma
Gmygamma = function(GFQzig1, GFQzig2){
  # GFQ of gamma
  Ggamma = sqrt(1/((1-GFQzig1$Gdelta)*GFQzig1$Ga))-sqrt(1/((1-GFQzig2$Gdelta)*GFQzig2$Ga))
  return(Ggamma)
}

ci.f = function(zigamma1,zigamma2,nr,alpha){
  #GFQ of gamma (=CV1-CV2)
  GFQzig1=GFQzig(zigamma1$ni0,zigamma1$ni1,zigamma1$zigdata1,nr=nr)
  GFQzig2=GFQzig(zigamma2$ni0,zigamma2$ni1,zigamma2$zigdata1,nr=nr)
  Ggamma=Gmygamma(GFQzig1, GFQzig2)
  fci=sort(Ggamma)
  l=fci[nr*alpha/2];u=fci[nr*(1-alpha/2)]
  return(list(lower=l,upper=u))
}

