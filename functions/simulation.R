source('functions/zigamma_data.R')

#--------simu----------
# simu contains all four methods
simu = function(n1, delta1, a1, b1, n2, delta2, a2, b2, trial,alpha){
  
  p_pb = 0; len_pb = 0; le_pb = 0; ue_pb = 0
  p_f = 0; len_f = 0; le_f = 0; ue_f = 0
  p_m1 = 0; len_m1 = 0; le_m1 = 0; ue_m1 = 0
  p_m2 = 0; len_m2 = 0; le_m2 = 0; ue_m2 = 0
  
  CV = sqrt(1/((1-delta1)*a1))-sqrt(1/((1-delta2)*a2))
  
  for(i in 1:trial){
    
    zigamma1 = rzigamma(n = n1, delta = delta1, shape = a1, scale = b1)
    zigamma2 = rzigamma(n = n2, delta = delta2, shape = a2, scale = b2)
    
    # PB
    dci_pb = ci.pb(zigamma1, zigamma2, nr=1000, alpha)
    l_pb = dci_pb$lower
    u_pb = dci_pb$upper  
    
    len_pb = len_pb+(u_pb-l_pb)
    if(l_pb <= CV & CV <= u_pb) {p_pb = p_pb+1} # coverage_probability
    if(CV < l_pb) {le_pb = le_pb+1} # lower_error_probability 
    if(CV > u_pb) {ue_pb = ue_pb+1} # upper_error_probability
    
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
  
  # pb
  cp_pb = p_pb/trial
  al_pb = len_pb/trial # average length
  l_pb = le_pb/trial
  u_pb = ue_pb/trial
  
  # fiducial
  cp_f = p_f/trial
  al_f = len_f/trial # average length
  l_f = le_f/trial
  u_f = ue_f/trial
  
  # mover1
  cp_m1 = p_m1/trial
  al_m1 = len_m1/trial # average length
  l_m1 = le_m1/trial
  u_m1 = ue_m1/trial
  
  # mover2
  cp_m2 = p_m2/trial
  al_m2 = len_m2/trial # average length
  l_m2 = le_m2/trial
  u_m2 = ue_m2/trial

  
  return(c(n1 = n1, n2 = n2, delta1 = delta1, delta2 = delta2, a1 = a1, a2 = a2,
           cp_pb = cp_pb,
           al_pb = al_pb,
           l_pb = l_pb,
           u_pb = u_pb,
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



a_col = c(1,2,5,10,20) # the collection of a
d_col= c(0.1, 0.3, 0.6)
delta1 = d_col[2]; delta2 = delta1
n1 = 132; n2 = 132
# 1 ==============================================
si_1 = c()
for(i in a_col){
  si1 = simu(n1, delta1, a1 = 1, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             trial = 1000, alpha = 0.05) 
  si_1 = rbind(si_1, si1)
  print(i)
}


write.csv(si_1, file = " r1.csv",
          na = "NA",row.names = F,
          fileEncoding = "")
# 2 ==============================================
si_2 = c()
for(i in a_col[-seq(1)]){
  si1 = simu(n1, delta1, a1 = 2, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             trial = 1000, alpha = 0.05) 
  si_2 = rbind(si_2, si1)
  print(i)
}


write.csv(si_2, file = " r2.csv",
          na = "NA",row.names = F,
          fileEncoding = "")
# 3 ==============================================

si_3 = c()
for(i in a_col[-seq(1,2)]){
  si1 = simu(n1, delta1, a1 = 5, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             trial = 1000, alpha = 0.05) 
  si_3 = rbind(si_3, si1)
}


write.csv(si_3, file = " r3.csv",
          na = "NA",row.names = F,
          fileEncoding = "")
# 4 ==============================================

si_4 = c()
for(i in a_col[-seq(1,3)]){
  si1 = simu(n1, delta1, a1 = 10, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             trial = 1000, alpha = 0.05) 
  si_4 = rbind(si_4, si1)
}


write.csv(si_4, file = " r4.csv",
          na = "NA",row.names = F,
          fileEncoding = "")
# 5 ==============================================

si_5 = c()
si1 = simu(n1, delta1, a1 = 20, b1 = 1, 
           n2, delta2, a2 = 20, b2 = 1, 
           trial = 1000, alpha = 0.05) 
si_5 = rbind(si_5, si1)

write.csv(si_5, file = " r5.csv",
          na = "NA",row.names = F,
          fileEncoding = "")

