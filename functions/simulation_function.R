source('functions/zigamma_data.R')
source('functions/1_bootstrap_method.R')
source('functions/2_fiducial_method.R')
source('functions/3_MOVER_method_ver1.R')
source('functions/4_MOVER_method_ver2.R')

#--------simu----------
simu = function(n1, delta1, a1, b1, n2, delta2, a2, b2, n_trial, alpha){
  # simu uses all four methods to obtain confidence intervals first, and then calculate 
  # coverage probabilities, average lengths, lower error probabilities, and upper error probabilities.
  
  # Args:
  # n1 (sample size), delta1, a1, b1 are parameters of the first ZI-gamma distribution
  # n2 (sample size), delta2, a2, b2 are parameters of the second ZI-gamma distribution
  # n_trial: the number of trials
  # alpha: significance level
  
  # Returns:
  # information of this trial (n1, delta1, a1, b1, n2, delta2, a2, b2)
  # coverage probability, average length, lower error probability, and upper error probability of these methods
  
  p_bootstrap = 0; len_bootstrap = 0; le_bootstrap = 0; ue_bootstrap = 0
  p_f = 0; len_f = 0; le_f = 0; ue_f = 0
  p_m1 = 0; len_m1 = 0; le_m1 = 0; ue_m1 = 0
  p_m2 = 0; len_m2 = 0; le_m2 = 0; ue_m2 = 0
  
  CV = sqrt(1/((1-delta1)*a1))-sqrt(1/((1-delta2)*a2))
  
  for(i in 1:n_trial){
    # sample from two ZI-gamma distributions
    zigamma1 = rzigamma(n = n1, delta = delta1, shape = a1, scale = b1)
    zigamma2 = rzigamma(n = n2, delta = delta2, shape = a2, scale = b2)
    
    # bootstrap
    dci_bootstrap = ci.pb(zigamma1, zigamma2, nr=1000, alpha)
    l_bootstrap = dci_bootstrap$lower
    u_bootstrap = dci_bootstrap$upper  
    
    len_bootstrap = len_bootstrap+(u_bootstrap-l_bootstrap)
    if(l_bootstrap <= CV & CV <= u_bootstrap) {p_bootstrap = p_bootstrap+1} # coverage_probability
    if(CV < l_bootstrap) {le_bootstrap = le_bootstrap+1} # lower_error_probability 
    if(CV > u_bootstrap) {ue_bootstrap = ue_bootstrap+1} # upper_error_probability
    
    # Fiducial
    dci_f=ci.f(zigamma1,zigamma2,nr=10000,alpha=alpha)
    l_f = dci_f$lower
    u_f = dci_f$upper
    
    len_f = len_f+(u_f-l_f)
    if(l_f <= CV & CV <= u_f) {p_f = p_f+1} # coverage_probability
    if(CV < l_f) {le_f = le_f+1} # lower_error_probability
    if(CV > u_f) {ue_f = ue_f+1} # upper_error_probability
    
    # mover1
    dci_m1 = ci.mover1(zigamma1, zigamma2, alpha)
    l_m1 = dci_m1$lower
    u_m1 = dci_m1$upper
    
    len_m1 = len_m1+(u_m1-l_m1)
    if(l_m1 <= CV & CV <= u_m1) {p_m1 = p_m1+1} # coverage_probability
    if(CV < l_m1) {le_m1 = le_m1+1} # lower_error_probability 
    if(CV > u_m1) {ue_m1 = ue_m1+1} # upper_error_probability

    # mover2
    dci_m2 = ci.mover2(zigamma1, zigamma2, alpha)
    l_m2 = dci_m2$lower
    u_m2 = dci_m2$upper  # after the confidence interval of delta obtained, perform log transformation to get CI of log(1-delta).dai

    len_m2 = len_m2+(u_m2-l_m2)
    if(l_m2 <= CV & CV <= u_m2) {p_m2 = p_m2+1} # coverage_probability
    if(CV < l_m2) {le_m2 = le_m2+1} # lower_error_probability
    if(CV > u_m2) {ue_m2 = ue_m2+1} # upper_error_probability
  }
  
  # bootstrap
  cp_bootstrap = p_bootstrap/n_trial # coverage probability
  al_bootstrap = len_bootstrap/n_trial # average length
  l_bootstrap = le_bootstrap/n_trial # lower error
  u_bootstrap = ue_bootstrap/n_trial # upper error
  
  # fiducial
  cp_f = p_f/n_trial
  al_f = len_f/n_trial # average length
  l_f = le_f/n_trial
  u_f = ue_f/n_trial
  
  # mover1
  cp_m1 = p_m1/n_trial
  al_m1 = len_m1/n_trial # average length
  l_m1 = le_m1/n_trial
  u_m1 = ue_m1/n_trial
  
  # mover2
  cp_m2 = p_m2/n_trial
  al_m2 = len_m2/n_trial # average length
  l_m2 = le_m2/n_trial
  u_m2 = ue_m2/n_trial

  
  return(c(n1 = n1, n2 = n2, delta1 = delta1, delta2 = delta2, a1 = a1, a2 = a2,
           cp_bootstrap = cp_bootstrap,
           al_bootstrap = al_bootstrap,
           l_bootstrap = l_bootstrap,
           u_bootstrap = u_bootstrap,
           cp_f = cp_f,
           al_f = al_f,
           l_f = l_f,
           u_f = u_f,
           cp_m1 = cp_m1,
           al_m1 = al_m1,
           l_m1 = l_m1,
           u_m1 = u_m1,
           cp_m2 = cp_m2,
           al_m2 = al_m2,
           l_m2 = l_m2,
           u_m2 = u_m2
           ))

}