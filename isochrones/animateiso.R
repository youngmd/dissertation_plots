# iso <- read.table('isoGC_BlueRed.dat', header=TRUE)
# 
# isoB <- subset(iso, Z == 0.000456)
# isoB$Group <- "Blue"
# 
# isoR <- subset(iso, Z == 0.00456)
# isoR$Group <- "Red"
# 
# iso <- rbind(isoB, isoR)
# intmag <- read.table('ssp03.dat', header=TRUE)
# intmag$age <- signif(log10(intmag$Age), 3)
# iso$age <- signif(iso$age, 3)
# intmag$BR <- signif(intmag$B - intmag$R, 3)
# 
# intmag <- data.frame(age=intmag$age, BR=intmag$BR)
# 
# iso2 <- merge(x=iso,y=intmag,by="age")
# iso2 <- subset(iso2, !is.na(BR))
library(ggplot2)
library(gganimate)

iso <- read.table('iso03.dat', header=TRUE)

isoB <- subset(iso, Z == 0.000456)
isoB$Group <- "Blue"

isoR <- subset(iso, Z == 0.00456)
isoR$Group <- "Red"

iso <- rbind(isoB, isoR)

p<-ggplot(iso, aes(x=B-R, y=V, frame=age))+
  geom_point(aes(color=Group), alpha=0.5) +
  scale_x_continuous(limits=c(-1, 5)) + 
  scale_y_reverse(limits=c(20, -15)) +
  theme_bw() +
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=12),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=10), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=10)) +
  scale_color_manual(values=c("#00B0F6",
                              "#F8766D"))

gg_animate(p, "isochrone.gif", interval=.1)

sspB <- read.table('ssp003.dat', header=TRUE)

sspR <- read.table('ssp03.dat', header=TRUE)

sspV <- read.table('ssp015.dat', header=TRUE)

ages <- data.frame(age=levels(factor(iso$age)),x=1,y=1)
sspB$Group <- "Blue"
sspR$Group <- "Red"
sspV$Group <- "Green"

ssp <- rbind(sspB, sspR, sspV)
ssp$age <- signif(log10(ssp$Age), 3)

p1 <- ggplot(data=ssp, aes(x=age,y=B-R, color=Group)) +
  geom_vline(aes(xintercept=age, frame=age), lty=2, color="black") +
  geom_line() +
  theme_bw() +
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=12),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=10), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=10)) +
  scale_color_manual(values=c("#00B0F6",
                              "#F8766D",
                              "green4"))

gg_animate(p1, "isometal.gif", interval=.1)

gcurves <- data.frame( BmR = 0, age=0, B=0, R=0)
sspB$age <- signif(log10(sspB$Age), 3)
sspR$age <- signif(log10(sspR$Age), 3)
for (i in 1:length(sspB$age)){
  B_BR = sspB[i,'B'] - sspB[i,'R']
  R_BR = sspR[i,'B'] - sspR[i,'R']
  age = sspB[i,'age'][1]
  dist <- data.frame(BmR = seq(from = -1, to = 1.8, length.out = 500), age=age)
  dist$B <- dnorm(dist$BmR,mean=B_BR, sd=0.1) * 0.7
  dist$R <- dnorm(dist$BmR,mean=R_BR, sd=0.1) * 0.3
  gcurves = rbind(gcurves, dist)
}

gcurves <- subset(gcurves, age != 0)

p2 <- ggplot(data=gcurves, aes(x=BmR, frame=age)) + 
  geom_ribbon(aes(ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  geom_ribbon(aes(ymin=0, ymax=R), color='black', fill='#F8766D', alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=12),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=10), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=10)) +
  ylab(expression('Density')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 

gg_animate(p2, "isocurves.gif", interval=.1)
