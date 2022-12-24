# simulation time
library("cowplot")
library("readxl")
library("ggplot2")
source('functions/zigamma_data.R')

f.simu.time<-function(n1, delta1, a1, b1, n2, delta2, a2, b2, trial, alpha){
  # Start the clock
  start_time = Sys.time()
  
  for(i in 1:trial){
    zigamma1=rzigamma(n=n1,delta=delta1,shape=a1,scale=b1)
    zigamma2=rzigamma(n=n2,delta=delta2,shape=a2,scale=b2)
    myfci=ci.f(zigamma1,zigamma2,nr=10000,alpha=alpha)
    
  }
  # Stop the clock
  end_time = Sys.time()
  time = end_time-start_time
  return(c(time = time))     
}
m1.simu.time<-function(n1, delta1, a1, b1, n2, delta2, a2, b2, trial, alpha){
  # Start the clock
  start_time = Sys.time()
  
  for(i in 1:trial){
    zigamma1 = rzigamma(n = n1, delta = delta1, shape = a1, scale = b1)
    zigamma2 = rzigamma(n = n2, delta = delta2, shape = a2, scale = b2)
    
    dci = ci.mover1(zigamma1, zigamma2, alpha)
  }
  # Stop the clock
  end_time = Sys.time()
  time = end_time-start_time
  return(c(time = time))     
}
m2.simu.time<-function(n1, delta1, a1, b1, n2, delta2, a2, b2, trial, alpha){
  # Start the clock
  start_time = Sys.time()
  
  for(i in 1:trial){
    zigamma1 = rzigamma(n = n1, delta = delta1, shape = a1, scale = b1)
    zigamma2 = rzigamma(n = n2, delta = delta2, shape = a2, scale = b2)
    
    dci = ci.mover2(zigamma1, zigamma2, alpha)
  }
  # Stop the clock
  end_time = Sys.time()
  time = end_time-start_time
  return(c(time = time))     
}
pb.simu.time<-function(n1, delta1, a1, b1, n2, delta2, a2, b2, trial, alpha){
  # Start the clock
  start_time = Sys.time()
  
  for(i in 1:trial){
    zigamma1 = rzigamma(n = n1, delta = delta1, shape = a1, scale = b1)
    zigamma2 = rzigamma(n = n2, delta = delta2, shape = a2, scale = b2)
    dci = ci.pb(zigamma1, zigamma2, nr = 1000, alpha)
  }
  # Stop the clock
  end_time = Sys.time()
  time = end_time-start_time
  return(c(time = time))     
}

trial = c(1000, 3000, 6000, 9000)
# time1 n = 25
simu.time = function(){
  time.f = c()
  # time.m1 = c()
  # time.m2 = c()
  # time.pb = c()
  for (i in trial){
    t = pb.simu.time(n1 = 25, delta1 = 0.1, a1 = 1, b1 = 1, 
                    n2 = 25, delta2 = 0.1, a2 = 1, b2 = 1, 
                    trial = i, alpha = 0.05)
    time.f = c(time.f, t)
  }
  return(time.f)
}

t1 = simu.time()
t2 = simu.time()
t3 = simu.time()
t4 = simu.time()



# time2 n = 50
simu.time.f = function(){
  time = c()
  for (i in trial){
    t = f.simu.time(n1 = 50, delta1 = 0.1, a1 = 1, b1 = 1, 
                    n2 = 50, delta2 = 0.1, a2 = 1, b2 = 1, 
                    trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t50.1 = simu.time.f()
simu.time.m1 = function(){
  time = c()
  for (i in trial){
    t = m1.simu.time(n1 = 50, delta1 = 0.1, a1 = 1, b1 = 1, 
                     n2 = 50, delta2 = 0.1, a2 = 1, b2 = 1, 
                     trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t50.2 = simu.time.m1()
simu.time.m2 = function(){
  time = c()
  for (i in trial){
    t = m2.simu.time(n1 = 50, delta1 = 0.1, a1 = 1, b1 = 1, 
                     n2 = 50, delta2 = 0.1, a2 = 1, b2 = 1, 
                     trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t50.3 = simu.time.m2()
simu.time.pb = function(){
  time = c()
  for (i in trial){
    t = pb.simu.time(n1 = 50, delta1 = 0.1, a1 = 1, b1 = 1, 
                     n2 = 50, delta2 = 0.1, a2 = 1, b2 = 1, 
                     trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t50.4 = simu.time.pb()

data50 = data.frame(t50 = c(t50.1, t50.2, t50.3, t50.4))




# time 3 n = 132
simu.time.f = function(){
  time = c()
  for (i in trial){
    t = f.simu.time(n1 = 132, delta1 = 0.1, a1 = 1, b1 = 1, 
                     n2 = 132, delta2 = 0.1, a2 = 1, b2 = 1, 
                     trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t132.1 = simu.time.f()
simu.time.m1 = function(){
  time = c()
  for (i in trial){
    t = m1.simu.time(n1 = 132, delta1 = 0.1, a1 = 1, b1 = 1, 
                    n2 = 132, delta2 = 0.1, a2 = 1, b2 = 1, 
                    trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t132.2 = simu.time.m1()
simu.time.m2 = function(){
  time = c()
  for (i in trial){
    t = m2.simu.time(n1 = 132, delta1 = 0.1, a1 = 1, b1 = 1, 
                    n2 = 132, delta2 = 0.1, a2 = 1, b2 = 1, 
                    trial = i, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t132.3 = simu.time.m2()
simu.time.pb = function(){
  time = c()
  for (i in trial){
    t = pb.simu.time(n1 = 132, delta1 = 0.1, a1 = 1, b1 = 1, 
                    n2 = 132, delta2 = 0.1, a2 = 1, b2 = 1, 
                    trial = 1000, alpha = 0.05)
    time = c(time, t)
  }
  return(time)
}
t132.4 = simu.time.pb()

data132 = data.frame(t132 = c(t132.1, t132.2, t132.3, t132.4))


# plot

data.plot3 = read_excel("data/simulation_results/simulation_results.xlsx", 
                        sheet = "plot3")
p1 = ggplot(data.plot3[1:16,], aes(x = Trial, y = Time1, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = data.plot3[1:16,]$Trial)+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("The number of trial") + ylab("Time (s)") + # Set axis labels
  ggtitle("Time comparsion for CI with n=25") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
p2 = ggplot(data.plot3[17:32,], aes(x = Trial, y = Time1, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = data.plot3[17:32,]$Trial)+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("The number of trial") + ylab("Time (s)") + # Set axis labels
  ggtitle("Time comparsion for CI with n=50") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
p3 = ggplot(data.plot3[33:48,], aes(x = Trial, y = Time1, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = data.plot3[33:48,]$Trial)+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("The number of trial") + ylab("Time (s)") + # Set axis labels
  ggtitle("Time comparsion for CI with n=132") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
png('plots/time_consumption.png', width = 800, height = 400)
plot_grid(p1, p2, p3,
          ncol = 3, nrow = 1)
dev.off()
