library(ggplot2)
library(tidyr)
library(mclust)
library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)

gcc90<- read.csv('gcc90.2.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

d <- subset(gcc90, galaxy_name == 'NGC4406')
d$bmr <- d$BR

x <- d$bmr
binw <- 1.2 * IQR(x) / length(x)^(1/3)

# binw <- 0.05
dist <- data.frame( BmR = seq(from = min(x)*0.95, to = max(x)*1.05, length.out = 500))

m <- Mclust(d$bmr, G=1, modelName='E')
sd = sqrt(m$parameters$variance$sigmasq)

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd) * (m$parameters$pro[1] * length(x) * binw)

#dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd2) * (m$parameters$pro[2] * length(x) * binw)
#dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sd3) * (m$parameters$pro[3] * length(x) * binw)

g0 <- ggplot() + 
  geom_histogram(data=d, aes(x=bmr), fill="grey", binwidth=binw) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  ylab(expression('Number')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 

g1 <- ggplot() + 
  geom_histogram(data=d, aes(x=bmr), fill="grey", binwidth=binw) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  ylab(expression('Number')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 


m <- Mclust(d$bmr, G=2, modelName='E')
sd = sqrt(m$parameters$variance$sigmasq)

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd) * (m$parameters$pro[1] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd) * (m$parameters$pro[2] * length(x) * binw)

g2 <- ggplot() + 
  geom_histogram(data=d, aes(x=bmr), fill="grey", binwidth=binw) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=R),  color='black', fill='#FF67A4', alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  ylab(expression('Number')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 

m <- Mclust(d$bmr, G=3, modelName='E')
sd = sqrt(m$parameters$variance$sigmasq)

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd) * (m$parameters$pro[1] * length(x) * binw)
dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd) * (m$parameters$pro[2] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sd) * (m$parameters$pro[3] * length(x) * binw)

g3 <- ggplot() + 
  geom_histogram(data=d, aes(x=bmr), fill="grey", binwidth=binw) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=V),  color='black', fill='#00BA38', alpha=0.3) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=R),  color='black', fill='#FF67A4', alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  ylab(expression('Number')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 

m <- Mclust(d$bmr, G=4, modelName='V')
sd1 = sqrt(m$parameters$variance$sigmasq[1])
sd2 = sqrt(m$parameters$variance$sigmasq[2])
sd3 = sqrt(m$parameters$variance$sigmasq[3])
sd4 = sqrt(m$parameters$variance$sigmasq[4])

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd1) * (m$parameters$pro[1] * length(x) * binw)
dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd2) * (m$parameters$pro[2] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sd3) * (m$parameters$pro[3] * length(x) * binw)
dist$P <- dnorm(dist$BmR,mean=m$parameters$mean[4], sd=sd4) * (m$parameters$pro[4] * length(x) * binw)



g4 <- ggplot() + 
  geom_histogram(data=d, aes(x=bmr), fill="grey", binwidth=binw) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=V),  color='black', fill='#00BA38', alpha=0.3) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=R),  color='black', fill='#FF67A4', alpha=0.3) +
  geom_ribbon(data = dist, aes( BmR, ymin=0, ymax=P),  color='black', fill='#B79F00', alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  ylab(expression('Number')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 

ggsave('mmExamples.g0.pdf', g0, width=4, height=4)
ggsave('mmExamples.g1.pdf', g1, width=4, height=4)
ggsave('mmExamples.g2.pdf', g2, width=4, height=4)
ggsave('mmExamples.g3.pdf', g3, width=4, height=4)
ggsave('mmExamples.g4.pdf', g4, width=4, height=4)