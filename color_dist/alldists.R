library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(gtable)
library(plyr)

gcc90<- read.csv('gcc90.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

for(i in levels(gcc90$galaxy_name)) { 
  gccdata <- subset(gcc90, galaxy_name == i)
  n <- nrow(gccdata)
  if(n < 40){
    gcc90 <- subset(gcc90, galaxy_name != i)
  }
}
gcc90 <- arrange(gcc90, galaxy_name)

g <- ggplot(data=gcc90, aes(x=BR)) +
  geom_histogram(fill="grey", binwidth=1/20) +
  facet_wrap(~ galaxy_name, scales="free_y",ncol=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
ylab(expression('Number')) +
xlab(expression('B - R'))

ggsave('alldists.pdf', g, width=8, height=8)