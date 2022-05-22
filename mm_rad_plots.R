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
library(plyr)

# mwgc <-read.table('harris_clean.txt',header=TRUE)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc <- subset(mwgc, V.R != 0)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc$BR <- mwgc$B.V + mwgc$V.R - (mwgc$A_B - mwgc$A_R)
# mwgc <- subset(mwgc, BR < 1.6)

gcc90<- read.csv('gcc90.2.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

extent = 19
gcc90<- subset(gcc90, galcen_radius < extent)

d <- subset(gcc90, galaxy_name == 'NGC4406')

m <- Mclust(d$BR, modelNames='E', G=3)

d$class <- as.character(m$classification)

g <- ggplot(data=d, aes(x=bx,y=by)) + geom_point(aes(color=class, shape=class))

datapath <- "~/gcc/Files_for_GCS_Database/n4406/"
cnames <- c("r","N","sd","sderr","lsd","lsderrhi","lsderrlo","good_frac","r_d","NGC","Ncontam","contam_frac")
#cnames2 <- c("r","N","sd","sderr","lsd","lsderrhi","lsderrlo","good_frac","r_d")

rp = read.table(paste(datapath,'/radial_profile.dat',sep=''), col.names = cnames)

rmin = 0
newrp <- data.frame(r=0, class='A', sd=0) 
for(i in 1:(length(rp$r))){
  rmax = rp[i,'r'] + (rp[i+1,'r'] - rp[i,'r']) / 2.0
  d1 <- subset(d, galcen_radius > rmin & galcen_radius < rmax)
  rmin <- rmax
  if(length(d1$galcen_radius) == 0){
    next
  }
  freqs <- count(d1, 'class')
  freqs$frac <- freqs$freq / sum(freqs$freq)
  freqs$sd <- freqs$frac * rp[i,'sd']
  freqs$r <- rp[i,'r']
  newrpt <- data.frame(r=freqs$r, class=freqs$class, sd=freqs$sd)
  newrp <- rbind(newrp, newrpt)
}

newrp <- subset(newrp, class != 'A')

newrpb <- subset(newrp, class=='1')
newrpb$sdr <- newrpb$sd / max(newrpb$sd)
newrpg <- subset(newrp, class=='2')
newrpg$sdr <- newrpg$sd / max(newrpg$sd)
newrpr <- subset(newrp, class=='3')
newrpr$sdr <- newrpr$sd / max(newrpr$sd)

newrp <- rbind(newrpb, newrpg, newrpr)

g2 <- ggplot(data=newrp, aes(x=r, y=sd)) + 
  geom_point(aes(color=class, shape=class)) + 
  geom_line(aes(color=class)) +
  theme_bw() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("blue3", "green3",
                              "red3"))

g3 <- ggplot(data=newrp, aes(x=r/extent, y=sdr)) + 
  geom_point(aes(color=class, shape=class)) + 
  geom_line(aes(color=class)) +
  theme_bw() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("blue3", "green3",
                              "red3"))

x <- d$BR
binw <- 1.2 * IQR(x) / length(x)^(1/3)
# binw <- 0.05
dist <- data.frame( BmR = seq(from = min(x)*0.95, to = max(x)*1.05, length.out = 500))

sd = sqrt(m$parameters$variance$sigmasq)

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd) * (m$parameters$pro[1] * length(x) * binw)
dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd) * (m$parameters$pro[2] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sd) * (m$parameters$pro[3] * length(x) * binw)

g1 <- ggplot() + 
  geom_histogram(data=d, aes(x=BR), fill="grey", binwidth=binw) +
  geom_ribbon(data = dist, aes(x=BmR, ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  geom_ribbon(data = dist, aes(x=BmR, ymin=0, ymax=R), color='black', fill='#FF67A4', alpha=0.3) +
  geom_ribbon(data = dist, aes(x=BmR, ymin=0, ymax=V), color='black', fill='#00BA38', alpha=0.3) +
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


db <- subset(d, class=='1')
dg <- subset(d, class=='2')
dr <- subset(d, class=='3')

g4 <- ggplot(data=db, aes(x=bx, y=by)) + geom_point(color='blue4',shape=0)
g5 <- ggplot(data=dg, aes(x=bx, y=by)) + geom_point(color='green4',shape=1)
g6 <- ggplot(data=dr, aes(x=bx, y=by)) + geom_point(color='red4',shape=2)

gplot <- arrangeGrob(g1, g2, g3, g4, g5, g6, nrow=2, ncol=3)
ggsave('mm_rad.pdf', gplot, width=8, height=6)
