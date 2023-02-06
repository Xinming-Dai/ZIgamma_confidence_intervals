# This is the simulation of MOVER method.
# CIs of alpha were constructed by Fiducial method.
# The estimator of a is constructed by MLE. losed-form estimation from 【Gamma distribution - Wikipedia】.
# ci.mover1 is the confidence interval function

library(EnvStats)
source('functions/zigamma_data.R')

# CI of a, which was constructed by Fiducial method.
a.ci = function(zigamma, nr, alpha){
  ni1 = zigamma$ni1
  zigdata1 = zigamma$zigdata1
  # GFQ of a
  y=zigdata1^(1/3)
  u=mean(y)
  s2=var(y)
  Z = rnorm(nr, 0,1); k = (ni1-1)/rchisq(nr, df = ni1-1)
  Gmu=u+Z*sqrt(s2/ni1)*sqrt(k)
  Gsigma2=s2*k
  Ga=1/9*((1+.5*Gmu^2/Gsigma2)+((1+.5*Gmu^2/Gsigma2)^2-1)^(1/2))
  fci=sort(Ga)
  l=fci[nr*alpha/2];u=fci[nr*(1-alpha/2)]
  return(list(lower=l,upper=u))
}

delta_.ci = function(zigamma, alpha){
  n = zigamma$ni0+zigamma$ni1
  delta_h = zigamma$ni0/n; q0 = 1-delta_h
  z0 = qnorm(1-alpha/2)
  delta_l = (delta_h+(z0^2/n)/2-sqrt(delta_h*q0*(z0^2/n)+(z0^2/n)^2/4))/(1+z0^2/n) # Wilson (1927)
  delta_u = (delta_h+(z0^2/n)/2+sqrt(delta_h*q0*(z0^2/n)+(z0^2/n)^2/4))/(1+z0^2/n)
  l = log(1-delta_u)
  u = log(1-delta_l)
  return(list(lower = l, upper = u))
}

# the confidence interval of CV
cv.ci1 = function (zigamma, alpha){
  # 用 MOVER 求 logcv 的置信区间
  n = zigamma$ni0+zigamma$ni1; e_delta = zigamma$ni0/n; theta1 = log(1-e_delta)
  es = egamma(zigamma$zigdata1, method = "mle"); e_a = es$parameters[1]; theta2 = log(e_a)
  
  theta1.ci = delta_.ci(zigamma, alpha); l1 = theta1.ci$lower; u1 = theta1.ci$upper
  theta2.ci = a.ci(zigamma, 10000, alpha); l2 = log(theta2.ci$lower); u2 = log(theta2.ci$upper)
  
  # MOVER
  L1 = theta1+theta2-sqrt((theta1-l1)^2+(theta2-l2)^2)
  U1 = theta1+theta2+sqrt((u1-theta1)^2+(u2-theta2)^2)
  L = exp(-1/2*U1); U =exp(-1/2*L1) # 转换成 CV 的置信区间
  return(list(e_delta = e_delta, e_a = e_a, lower = L, upper = U ))
  
}

# the confidence interval of the difference of CVs
ci.mover1 = function(zigamma1, zigamma2, alpha){
  # utilize MOEVR to construct the confidence interval of the difference of CVs
  
  # Args:
  # zigamma1: 1st sample
  # zigamma2: 2nd sample
  # nr: number of resample size
  # alpha: significance level
  
  # Returns:
  # lower: the lower bound of the confidence interval
  # upper: the upper bound of the confidence interval
  
  ci1 = cv.ci1(zigamma1, alpha); l1 = ci1$lower; u1 = ci1$upper
  ci2 = cv.ci1(zigamma2, alpha); l2 = ci2$lower; u2 = ci2$upper
  
  theta1 = 1/sqrt((1-ci1$e_delta)*ci1$e_a)
  theta2 = 1/sqrt((1-ci2$e_delta)*ci2$e_a)
  
  
  L = theta1-theta2-sqrt((theta1-l1)^2+(u2-theta2)^2)
  U = theta1-theta2+sqrt((u1-theta1)^2+(theta2-l2)^2)
  names(L) = NULL; names(U) = NULL
  return(list(lower = L, upper = U))
}

