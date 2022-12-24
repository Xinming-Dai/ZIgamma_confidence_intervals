# This is the simulation of MOVER method.
# Ver2. CIs of scale parameters, a, were constructed by Score CI【3-2.1】.
# ci.mover2 is the confidence interval function

library(EnvStats) # MLE of gamma
source('functions/zigamma_data.R')

# alpha = 0.05; zigamma = zigamma1
# CI for CV of gammma which was constructed by Score CI.
CV2.ci = function(zigamma, alpha){
  ni1 = zigamma$ni1; x = zigamma$zigdata1
  z1 = ni1*log(mean(x))-sum(log(x))
  z = qnorm(1-alpha/2)
  # e_a = ni1*sum(x)/(ni1*sum(x*log(x))-sum(log(x))*sum(x)) # From Gamma Distribution Wikipedia
  e_a = 1/(2*(log(mean(x))-sum(log(x))/ni1))
  
  # Score CI
  l = sqrt(2/ni1*(z1-z*sqrt(ni1/(2*e_a^2))))
  u = sqrt(2/ni1*(z1+z*sqrt(ni1/(2*e_a^2))))
  return(list(lower = l, upper = u))
}

# delta_.ci is the CI of -1/2*log(1-delta)
delta_.ci = function(zigamma, alpha){
  n = zigamma$ni0+zigamma$ni1
  delta_h = zigamma$ni0/n; q0 = 1-delta_h
  z0 = qnorm(1-alpha/2)
  delta_l = (delta_h+(z0^2/n)/2-sqrt(delta_h*q0*(z0^2/n)+(z0^2/n)^2/4))/(1+z0^2/n) # Wilson (1927)
  delta_u = (delta_h+(z0^2/n)/2+sqrt(delta_h*q0*(z0^2/n)+(z0^2/n)^2/4))/(1+z0^2/n)
  l = -1/2*log(1-delta_l)
  u = -1/2*log(1-delta_u)
  return(list(lower = l, upper = u))
}

# the confidence interval of CV of ZIG
cv.ci2 = function (zigamma, alpha){
  # 用 MOVER 求 logcv 的置信区间
  
  # Binomial 
  n = zigamma$ni0+zigamma$ni1
  e_delta = zigamma$ni0/n
  theta1 = -1/2*log(1-e_delta)
  
  # Gamma
  es = egamma(zigamma$zigdata1, method = "mle")
  e_a = es$parameters[1]
  theta2 = -1/2*log(e_a)
  
  theta1.ci = delta_.ci(zigamma, alpha); l1 = theta1.ci$lower; u1 = theta1.ci$upper
  theta2.ci = CV2.ci(zigamma, alpha); l2 = log(theta2.ci$lower); u2 = log(theta2.ci$upper)
  
  # MOVER
  L1 = theta1+theta2-sqrt((theta1-l1)^2+(theta2-l2)^2)
  U1 = theta1+theta2+sqrt((u1-theta1)^2+(u2-theta2)^2)
  L = exp(L1); U = exp(U1) # 转换成 CV 的置信区间

  return(list(e_delta = e_delta, e_a = e_a, lower = L, upper = U))
  
}


# the confidence interval of the difference of CVs
ci.mover2 = function(zigamma1, zigamma2, alpha){
  # utilize MOEVR to construct the confidence interval of the difference of CVs
  
  ci1 = cv.ci2(zigamma1, alpha); l1 = ci1$lower; u1 = ci1$upper
  ci2 = cv.ci2(zigamma2, alpha); l2 = ci2$lower; u2 = ci2$upper
  
  theta1 = 1/sqrt((1-ci1$e_delta)*ci1$e_a)
  theta2 = 1/sqrt((1-ci2$e_delta)*ci2$e_a)
  
  L = theta1-theta2-sqrt((theta1-l1)^2+(u2-theta2)^2)
  U = theta1-theta2+sqrt((u1-theta1)^2+(theta2-l2)^2)
  names(L) = NULL; names(U) = NULL
  return(list(lower = L, upper = U))
}

