# comparing different methods
library("ggplot2")
library("readxl")
library("cowplot")

data = read_excel("data/simulation_results/simulation_results.xlsx", 
                  sheet = "plot1")
data$group = as.factor(data$group)
p1 = ggplot(data, aes(x = nratio, y = cp, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("Sample size") + ylab("Coverage probailities") + # Set axis labels
  ggtitle("Plot of coverage probailities by sample size") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
  
p2 = ggplot(data, aes(x = nratio, y = al, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("Sample size") + ylab("Average lengths") + # Set axis labels
  ggtitle("Plot of average lengths by sample size") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
#============================================================

data2 = read_excel("data/simulation_results/simulation_results.xlsx", 
                  sheet = "plot2")
p3= ggplot(data2, aes(x = dratio, y = cp, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("Probabilities of zero values") + ylab("Coverage probailities") + # Set axis labels
  ggtitle("Plot of coverage probailities by probabilities of zero values") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))
p4 = ggplot(data2, aes(x = dratio, y = al, group = group, colour=group))+
  geom_line()+
  geom_point()+
  scale_colour_hue(name = "Methods")+      # Set legend title
  xlab("Probabilities of zero values") + ylab("Average lengths") + # Set axis labels
  ggtitle("Plot of average lengths by probabilities of zero values") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey80", size = 0.2),
        panel.border = element_rect(color = "grey80", fill = NA))

png('plots/result_plot.png', width = 800, height = 600)
plot_grid(p1, p2, p3, p4,
          ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
dev.off()
