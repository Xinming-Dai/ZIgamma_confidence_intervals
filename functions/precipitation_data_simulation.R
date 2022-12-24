# to gamma distribution test
# The sample size is 132.
library("readxl")
library("goft")
library("ggplot2")
source('functions/1_bootstrap_method.r')
source('functions/2_fiducial_method.R')
source('functions/3_MOVER_method_ver1.R')
source('functions/4_MOVER_method_ver2.R')
source('functions/zigamma_data.R')

# Beijing
B = read_excel("data/precipitation.xlsx", 
               sheet = "Beijing")
B = B[, -1]
Beijing = c()
for(i in 1:dim(B)[1]){
  Beijing = c(Beijing, unlist(B[i,]))
}
names(Beijing) = NULL

# Density plot of data
B = data.frame(Beijing = Beijing)
ggplot(B,aes(x=Beijing))+geom_density()
Beijing = zig(Beijing)

# estimation
delta1 = Beijing$ni0/132
e1 = gamma_fit(Beijing$zigdata1)

hb = hist(Beijing$zigdata)
ks.test(Beijing$zigdata1, "pgamma", 0.67)

gamma_test(Beijing$zigdata1)
# p-value = 0.9726

# it seems that the lognormal model is more suitable.
fit1 = glm(hb$counts~hb$breaks[-8], family = Gamma)
fit2 = glm(hb$counts~hb$breaks[-8], family=gaussian(link="log"))
AIC(fit1)
AIC(fit2)

# plot data and gamma
gamma = rgamma(1000, shape = e1['shape', ], scale = e1['scale', ])
gammaData = data.frame(gamma = gamma)
ggplot(gammaData, aes(x = gamma))+
  geom_histogram(aes(y =..density..),
                 colour = "black",
                 fill = "white") +
  stat_function(fun = dgamma, bins = 10, args = list(shape = 0.6787571, scale = 83.7041872))

#Zhengzhou===========================================
Z = read_excel("data/precipitation.xlsx", 
               sheet = "Zhengzhou")
Z = Z[, -1]

Zhengzhou = c()
for(i in 1:dim(Z)[1]){
  Zhengzhou = c(Zhengzhou, unlist(Z[i,]))
}
names(Zhengzhou) = NULL

# Density plot of data
myData = data.frame(Precipitation = c(B$Beijing, Zhengzhou), 
                    Location = as.factor(c(rep("Beijing", length(B$Beijing)),
                                           rep("Zhengzhou", length(Zhengzhou)))))
png('plots/precipitation_density.png', width = 600, height = 400)
ggplot(myData,aes(x = Precipitation, group = Location, colour = Location))+
  geom_density()+
  scale_colour_hue(name = "Methods")+
  ggtitle("Density plot of Beijing and Zhengzhou precipitation data")+
  scale_colour_brewer(palette = "Set1")+
  ylab("Density")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
dev.off()

# Zig(Zhengzhou)
Zhengzhou = zig(Zhengzhou)

# estimation
delta2 = Zhengzhou$ni0/132
e2 = gamma_fit(Zhengzhou$zigdata1)

hz = hist(Zhengzhou$zigdata)
gamma_test(Zhengzhou$zigdata)
# p-value = 0.4429

# it seems that the lognormal model is more suitable.
fit1 = glm(hz$counts~hz$breaks[-7], family = Gamma)
fit2 = glm(hz$counts~hz$breaks[-7], family=gaussian(link="log"))
AIC(fit1)
AIC(fit2)
# 4 method test===============
sqrt(1/((1-delta1)*e1[1]))-sqrt(1/((1-delta2)*e2[1]))
zigamma1 = Beijing
zigamma2 = Zhengzhou
alpha = 0.05
dci1 = ci.pb(zigamma1, zigamma2, nr = 1000, alpha)
dci2 = ci.f(zigamma1,zigamma2,nr=10000,alpha=alpha)
dci3 = ci.mover1(zigamma1, zigamma2, alpha)
dci4 = ci.mover2(zigamma1, zigamma2, alpha)


