# this file uses `simu` function from simulation_function.R file to simulate
source('functions/simulation_function.R')
set.seed(2023)

a_col = c(1,2,5,10,20) # the collection of a
d_col= c(0.1, 0.3, 0.6) # the collection of b
delta1 = d_col[2]; delta2 = delta1
n1 = 132; n2 = 132
# 1 ==============================================
simulation_result_1 = c()
for(i in a_col){
  si1 = simu(n1, delta1, a1 = 1, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             n_trial = 1000, alpha = 0.05) 
  simulation_result_1 = rbind(simulation_result_1, si1)
  print(paste("finished", i, "th a."))
}

write.csv(simulation_result_1, file = " r1.csv",
          na = "NA",row.names = F,
          fileEncoding = "")
# 2 ==============================================
simulation_result_2 = c()
for(i in a_col[-seq(1)]){
  si1 = simu(n1, delta1, a1 = 2, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             n_trial = 1000, alpha = 0.05) 
  simulation_result_2 = rbind(simulation_result_2, si1)
  print("finished", i, "th a.")
}

write.csv(simulation_result_2, file = " r2.csv",
          na = "NA",row.names = F,
          fileEncoding = "")

# 3 ==============================================
simulation_result_3 = c()
for(i in a_col[-seq(1,2)]){
  si1 = simu(n1, delta1, a1 = 5, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             n_trial = 1000, alpha = 0.05) 
  simulation_result_3 = rbind(simulation_result_3, si1)
  print("finished", i, "th a.")
}

write.csv(simulation_result_3, file = " r3.csv",
          na = "NA",row.names = F,
          fileEncoding = "")

# 4 ==============================================
simulation_result_4 = c()
for(i in a_col[-seq(1,3)]){
  si1 = simu(n1, delta1, a1 = 10, b1 = 1, 
             n2, delta2, a2 = i, b2 = 1, 
             n_trial = 1000, alpha = 0.05) 
  simulation_result_4 = rbind(simulation_result_4, si1)
  print("finished", i, "th a.")
}

write.csv(simulation_result_4, file = " r4.csv",
          na = "NA",row.names = F,
          fileEncoding = "")

# 5 ==============================================
simulation_result_5 = c()
si1 = simu(n1, delta1, a1 = 20, b1 = 1, 
           n2, delta2, a2 = 20, b2 = 1, 
           n_trial = 1000, alpha = 0.05) 
simulation_result_5 = rbind(simulation_result_5, si1)

write.csv(simulation_result_5, file = " r5.csv",
          na = "NA",row.names = F,
          fileEncoding = "")

