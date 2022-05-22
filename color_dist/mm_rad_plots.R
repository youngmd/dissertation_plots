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
library(MASS)
library(cluster)
library(flux)


get_hd_ellipse <- function(xy, coverage) {
  fit <- CovMve(xy, alpha = coverage)
  points_in_ellipse <- xy[fit@best, ]
  ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
  ellipse_boundary
}

get_total <- function(df){
     binw <- df[2,'x'] - df[1,'x']
     df$area <- binw * 10^df$y
     df$total <- cumsum(df$area)
     reff <- approx(df$total, df$x, xout=0.5*max(df$total))
     reff$y
}

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

azimuth <- function(x, y)
{
  t <- atan(y/x)
  t <- t * 180/pi
  t <- t + 90
  t
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 500){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)     
  return(data.frame(x = xx, y = yy)) 
}

ellipsefun <- function(center = c(0,0), pa = 0, a = 1, b= 0.5, npoints = 500){
  tt <- seq(0,2*pi,length.out = npoints)
  xx = center[1] + a * cos(tt) * cos(pa) - b * sin(tt) * sin(pa)
  yy = center[2] + b * sin(tt) * cos(pa) + a * cos(tt) * sin(pa)
  return(data.frame(x = xx, y = yy)) 
}

d$azimuth <- azimuth(d$bx - galcenx, d$by - galceny)
d2$azimuth <- azimuth(d2$bx - galcenx, d2$by - galceny)
# mwgc <-read.table('harris_clean.txt',header=TRUE)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc <- subset(mwgc, V.R != 0)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc$BR <- mwgc$B.V + mwgc$V.R - (mwgc$A_B - mwgc$A_R)
# mwgc <- subset(mwgc, BR < 1.6)

m <- Mclust(d$BR, modelNames='E', G=3)
errs <- MclustBootstrap(m, type="bs")
errs <- summary(errs, what="se")
lerrs <- errs$pro
meanl <- mean(lerrs)

mixerr <- mean(errs$pro)

d$class <- as.character(m$classification)
d$uncert <- m$uncertainty
#d <- subset(d, uncert < 0.49)
scale <- mean(sqrt((d$bx - galcenx)^2 + (d$by - galceny)^2) / d$galcen_radius)
extentcirc <- circleFun(c(galcenx, galceny), 2*extent*scale, npoints=100)
galellip <- ellipsefun(c(galcenx, galceny), pa=(galpa+90)*(pi/180), a=galmajor*scale, b=galminor*scale, npoints=200)

rmin = 0
newrp <- data.frame(r=0, class='A', sd=0, sderr=0, frac=0) 
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
  freqs$sderr <- sqrt((mixerr/freqs$frac)^2 + (rp[i,'sderr']/rp[i,'sd'])^2)*freqs$sd 
  freqs$r <- rp[i,'r']
  newrpt <- data.frame(r=freqs$r, class=freqs$class, sd=freqs$sd, sderr=freqs$sderr, frac=freqs$frac)
  newrp <- rbind(newrp, newrpt)
}


newrp$lsd <- log10(newrp$sd)
newrp$lsderrhi <- log10(newrp$sd + newrp$sderr) - newrp$lsd

newrp <- subset(newrp, class != 'A')

newrp$dr <- newrp$r^0.25

newrpb <- subset(newrp, class=='1')
newrpb$sdr <- newrpb$sd / max(newrpb$sd)
newrpg <- subset(newrp, class=='2')
newrpg$sdr <- newrpg$sd / max(newrpg$sd)
newrpr <- subset(newrp, class=='3')
newrpr$sdr <- newrpr$sd / max(newrpr$sd)

newrp <- rbind(newrpb, newrpg, newrpr)

pfb <- summary(lm(log10(sd) ~ log10(r), data=newrpb, weight=1.0 / lsderrhi^2))
dfb <- summary(lm(log10(sd) ~ dr, data=newrpb, weight=1.0 / lsderrhi^2))

print(pfb)
print(dfb)

pfb$a0 <- pfb$coefficients[1,1]
pfb$a1 <- pfb$coefficients[2,1]
dfb$a0 <- dfb$coefficients[1,1]
dfb$a1 <- dfb$coefficients[2,1]

pfv <- summary(lm(log10(sd) ~ log10(r), data=newrpg, weight=1.0 / lsderrhi^2))
dfv <- summary(lm(log10(sd) ~ dr, data=newrpg, weight=1.0 / lsderrhi^2))

print(pfv)
print(dfv)

pfv$a0 <- pfv$coefficients[1,1]
pfv$a1 <- pfv$coefficients[2,1]
dfv$a0 <- dfv$coefficients[1,1]
dfv$a1 <- dfv$coefficients[2,1]

pfr <- summary(lm(log10(sd) ~ log10(r), data=newrpr, weight=1.0 / lsderrhi^2))
dfr <- summary(lm(log10(sd) ~ dr, data=newrpr, weight=1.0 / lsderrhi^2))
print(pfr)
print(dfr)

pfr$a0 <- pfr$coefficients[1,1]
pfr$a1 <- pfr$coefficients[2,1]
dfr$a0 <- dfr$coefficients[1,1]
dfr$a1 <- dfr$coefficients[2,1]

fitlines <- data.frame( r = seq(from = min(newrp$r), 
                                to = max(newrp$r), 
                                length.out = 500))

fitlineb <- fitlines
fitlineb$pow <- pfb$a0 + pfb$a1 * log10(fitlineb$r)
fitlineb$dev <- dfb$a0 + dfb$a1 * fitlineb$r^(0.25)

fitlinev <- fitlines
fitlinev$pow <- pfv$a0 + pfv$a1 * log10(fitlinev$r)
fitlinev$dev <- dfv$a0 + dfv$a1 * fitlinev$r^(0.25)

fitliner <- fitlines
fitliner$pow <- pfr$a0 + pfr$a1 * log10(fitliner$r)
fitliner$dev <- dfr$a0 + dfr$a1 * fitliner$r^(0.25)

newrp[newrp$sderr > newrp$sd,'sderr'] <- newrp[newrp$sderr > newrp$sd,'sd']

fits <- data.frame(pop='a', bestfit = 'a', a0=0, a1=0, a0err=0, a1err=0, r10=0)

if(pfb$fstatistic[1] < dfb$fstatistic[1]){
  fitlineb$pow <- fitlineb$dev
  f <- data.frame(pop='B', bestfit = 'dev', a0=dfb$a0, a1=dfb$a1, a0err=dfb$coefficients[1,2], a1err=dfb$coefficients[1,2])
} else {
  f <- data.frame(pop='B', bestfit = 'pow', a0=pfb$a0, a1=pfb$a1, a0err=pfb$coefficients[1,2], a1err=pfb$coefficients[1,2])
  
}

t <- data.frame(x=fitlineb$r, y=fitlineb$pow)
f$r10 <- get_total(t)
fits <- rbind(fits, f)

if(pfv$fstatistic[1] < dfv$fstatistic[1]){
  fitlinev$pow <- fitlinev$dev
  f <- data.frame(pop='V', bestfit = 'dev', a0=dfv$a0, a1=dfv$a1, a0err=dfv$coefficients[1,2], a1err=dfv$coefficients[1,2])
} else {
  f <- data.frame(pop='V', bestfit = 'pow', a0=pfv$a0, a1=pfv$a1, a0err=pfv$coefficients[1,2], a1err=pfv$coefficients[1,2])
}

t <- data.frame(x=fitlinev$r, y=fitlinev$pow)
f$r10 <- get_total(t) 
fits <- rbind(fits, f)

if(pfr$fstatistic[1] < dfr$fstatistic[1]){
  fitliner$pow <- fitliner$dev
  f <- data.frame(pop='R', bestfit = 'dev', a0=dfr$a0, a1=dfr$a1, a0err=dfr$coefficients[1,2], a1err=dfr$coefficients[1,2])
} else {
  f <- data.frame(pop='R', bestfit = 'pow', a0=pfr$a0, a1=pfr$a1, a0err=pfr$coefficients[1,2], a1err=pfr$coefficients[1,2])
}

t <- data.frame(x=fitliner$r, y=fitliner$pow)
f$r10 <- get_total(t)

fits <- rbind(fits, f)
fits$r10 <- fits$r10
fits <- subset(fits, a0 != 0)

fits$galaxy <- galname

outer = max(newrp$r)
g2 <- ggplot(data=newrp, aes(x=r, y=sd)) + 
  geom_point(aes(color=class, shape=class), size=0.8) + 
  geom_errorbar(aes(color=class, ymin=sd-sderr,ymax=sd+sderr), width=0.02*log10(outer), size=0.2) +
  geom_path(data = fitlineb, aes( r, 10^pow), linetype=1, size=0.3, color='#00B0F6') +
  geom_path(data = fitlinev, aes( r, 10^pow), linetype=2, size=0.3, color='#00BA38') +
  geom_path(data = fitliner, aes( r, 10^pow), linetype=3, size=0.3, color='#F8766D') +
  #geom_smooth(aes(color=class), method="glm") +
  ylab(expression(Number/arcmin^2)) + 
  xlab(expression('Projected radial distance (arcmin)')) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=8),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=6), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=6)) +
  scale_color_manual(values=c("#00B0F6", "#00BA38",
                              "#F8766D")) +
  scale_y_continuous(trans='asinh',breaks=c(0,1,10,100),labels=c("0",expression(10^0),expression(10^1),expression(10^2)),limits=c(-0.5,NA)) +
  scale_x_log10(breaks=c(1,2,5,10,20,50,100))

g1 <- ggplot(data=newrp, aes(x=r, y=sd)) + 
  geom_point(aes(color=class, shape=class), size=0.8) + 
  geom_errorbar(aes(color=class, ymin=sd-sderr,ymax=sd+sderr), width=0.02*outer, size=0.2) +
  geom_path(data = fitlineb, aes( r, 10^pow), linetype=1, size=0.3, color='#00B0F6') +
  geom_path(data = fitlinev, aes( r, 10^pow), linetype=2, size=0.3, color='#00BA38') +
  geom_path(data = fitliner, aes( r, 10^pow), linetype=3, size=0.3, color='#F8766D') +
  #geom_smooth(aes(color=class), method="glm") +
  ylab(expression(Number/arcmin^2)) + 
  xlab(expression('Projected radial distance (arcmin)')) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=8),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=6), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=6)) +
  scale_color_manual(values=c("#00B0F6", "#00BA38",
                              "#F8766D")) +
  annotate("text", x=Inf,y=Inf
           , label = galname,size = 4, hjust = 1.1, vjust = 2.0)


db <- subset(d, class=='1')
d2b <- subset(d2, BR < max(db$BR))
d2b$class <- '1'

dr <- subset(d, class=='3')
d2r <- subset(d2, BR > min(dr$BR))
d2r$class <- '3'

dg <- subset(d, class=='2')
d2g <- subset(d2, BR < min(dr$BR) & BR > max(db$BR))
d2g$class <- '2'


#daz <- subset(d, uncert < 0.45)
daz <- rbind(d2b, d2g, d2r)
daz <- subset(daz, galcen_radius < max(galmajor,reff))

x <- daz$azimuth
binw <- 1.5 * IQR(x) / length(x)^(1/3)

bins <- round((max(daz$azimuth) - min(daz$azimuth)) / binw)
binw <- (max(daz$azimuth) - min(daz$azimuth)) / bins

azdist <- data.frame(az=999, class='A', freq=0, frac=0, fracall=0, freqerr=0, fracallerr=0)
for(i in 1:bins){
  minaz <- (i-1)*binw + min(d$azimuth)
  maxaz <- minaz + binw
  d1 <- subset(daz, azimuth < maxaz & azimuth >= minaz)
  freqs <- count(d1, 'class')
  freqs2 <- count(daz,'class')
  freqs2 <- subset(freqs2, (freqs2$class %in% freqs$class) )
  freqs$frac <- freqs$freq / sum(freqs$freq)
  freqs$freqerr <- sqrt(freqs$freq)
  freqs$fracerr <- mixerr
  freqs$fracall <- freqs$freq / freqs2$freq 
  freqs$fracallerr <- (sqrt(freqs$freq)/freqs$freq + sqrt(freqs2$freq)/freqs2$freq)*freqs$fracall
  azdist2 <- data.frame(az=minaz +0.5*binw, class=freqs$class, freq=freqs$freq, frac=freqs$frac, fracall=freqs$fracall, freqerr=freqs$freqerr, fracallerr=freqs$fracallerr)
  azdist <- rbind(azdist, azdist2)
}

azdist <- subset(azdist, class != 'A')
azdist2 <- azdist
azdist2$az <- azdist$az - 0.5*binw
azdist4 <- subset(azdist, az > 0.95*max(az))
azdist4$az <- azdist4$az + 0.5*binw

azdist3 <- rbind(azdist,azdist2)
azdist3 <- rbind(azdist3, azdist4)

g3a <- ggplot(data=azdist3, aes(x=az, y=fracall)) +
  geom_step(aes(color=class, linetype=class), size=0.3) + 
  geom_point(data=azdist, aes(x=az, y=fracall, color=class, shape=class), size=0.8) +
  geom_errorbar(data=azdist, aes(x=az, ymin=fracall-fracallerr, ymax=fracall+fracallerr, color=class), width=3, size=0.3) +
  geom_vline(aes(xintercept=galpa), color='grey20', size=0.3, linetype=4) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=8),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=6),
        axis.text.x=element_text(margin=margin(5,0,0,0), size=6)) +
  scale_color_manual(values=c("#00B0F6", "#00BA38","#F8766D")) +
  scale_linetype_manual(values=c(1,2,3)) +
  scale_y_continuous(limits=c(0.9*min(azdist$fracall-azdist$fracallerr),NA)) +
  scale_x_continuous(expand=c(0,0)) +
  xlab(expression(paste("Position angle (",degree,")"))) +
  ylab("Fraction")

g3b <- ggplot(data=azdist3, aes(x=az, y=freq)) +
  geom_step(aes(color=class, linetype=class), size=0.3) + 
  geom_point(data=azdist, aes(x=az, y=freq, color=class, shape=class), size=0.8) +
  geom_errorbar(data=azdist, aes(x=az, ymin=freq-freqerr, ymax=freq+freqerr, color=class), width=3, size=0.3) +
  geom_vline(aes(xintercept=galpa), color='grey20', size=0.3, linetype=4) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=8),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=6),
        axis.text.x=element_text(margin=margin(5,0,0,0), size=6)) +
  scale_color_manual(values=c("#00B0F6", "#00BA38","#F8766D")) +
  scale_linetype_manual(values=c(1,2,3)) +
  scale_y_continuous(limits=c(0.9*min(azdist$frac),NA)) +
  scale_x_continuous(expand=c(0,0)) +
  xlab(expression(paste("Position angle (",degree,")"))) +
  ylab(expression(N[GC]))

g3 <- arrangeGrob(g3a, g3b,  nrow=2, ncol=1)

# g3 <- ggplot(data=newrp, aes(x=r/extent, y=sdr)) + 
#   geom_point(aes(color=class, shape=class), size=1.2) + 
#   geom_line(aes(color=class), size=0.2) +
#   ylab(expression(sigma/sigma[center])) + 
#   xlab(expression('Projected radial distance / '~GCS[extent])) +
  # theme_bw() +
  # theme(legend.position="none") +
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       legend.title=element_blank(),
  #       legend.key = element_blank(),
  #       axis.line = element_line(colour = "black"),
  #       axis.ticks.length=unit(-0.08, "cm"),
  #       axis.title=element_text(size=8),
  #       axis.text.y=element_text(margin=margin(0,5,0,0), size=6),
  #       axis.text.x=element_text(margin=margin(5,0,0,0), size=6)) +
#   scale_color_manual(values=c("#00B0F6", "#00BA38",
#                               "#F8766D"))

x <- d$BR
binw <- 1.2 * IQR(x) / length(x)^(1/3)
# binw <- 0.05
dist <- data.frame( BmR = seq(from = min(x)*0.95, to = max(x)*1.05, length.out = 500))

sd = sqrt(m$parameters$variance$sigmasq)

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd) * (m$parameters$pro[1] * length(x) * binw)
dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd) * (m$parameters$pro[2] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sd) * (m$parameters$pro[3] * length(x) * binw)

g0 <- ggplot() + 
  geom_histogram(data=d, aes(x=BR), fill="grey", binwidth=binw) +
  geom_ribbon(data = dist, aes(x=BmR, ymin=0, ymax=B), color='black', fill='#00B0F6', alpha=0.3) +
  geom_ribbon(data = dist, aes(x=BmR, ymin=0, ymax=R), color='black', fill='#F8766D', alpha=0.3) +
  geom_ribbon(data = dist, aes(x=BmR, ymin=0, ymax=V), color='black', fill='#00BA38', alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=8),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=6), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=6)) +
  ylab(expression('Number')) + 
  xlab(expression('B - R')) + 
  scale_y_continuous(expand= c(0,0)) 


minx = 0.99*min(d2$bx,mask$x1, mask$x2)
maxx = 1.01*max(d2$bx,mask$x1, mask$x2)
miny = 0.99*min(d2$by,mask$y1, mask$y2)
maxy = 1.01*max(d2$by,mask$y1, mask$y2)

g4 <- ggplot() + 
  geom_rect(data=mask, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.3) +
  geom_point(data=db, aes(x=bx, y=by), color='#00B0F6',shape=16) + 
  geom_point(data=d2b, aes(x=bx, y=by), shape=1, color='#00B0F6') +
  geom_path(data=extentcirc, aes(x=x, y=y), color='black', linetype=2, size=0.25) +
  geom_path(data=galellip, aes(x=x, y=y), color='black', linetype=1, size=0.25) +
  geom_point(aes(x=galcenx, y=galceny), color='black', shape=5) +
  scale_x_continuous(limits=c(minx,maxx)) + 
  scale_y_continuous(limits=c(miny,maxy))

g5 <- ggplot() + 
  geom_rect(data=mask, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.3) +
  geom_point(data=dg, aes(x=bx, y=by), color='#00BA38',shape=17) +
  geom_point(data=d2g, aes(x=bx, y=by), color='#00BA38',shape=2) + 
  geom_path(data=extentcirc, aes(x=x, y=y), color='black', linetype=2, size=0.25) +
  geom_path(data=galellip, aes(x=x, y=y), color='black', linetype=1, size=0.25) +
  geom_point(aes(x=galcenx, y=galceny), color='black', shape=5) + 
  scale_x_continuous(limits=c(minx,maxx)) + 
  scale_y_continuous(limits=c(miny,maxy))

g6 <- ggplot() + 
  geom_rect(data=mask, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.3) +
  geom_point(data=dr, aes(x=bx, y=by), color='#F8766D',shape=15) + 
  geom_point(data=d2r, aes(x=bx, y=by), color='#F8766D',shape=0) + 
  geom_path(data=extentcirc, aes(x=x, y=y), color='black', linetype=2, size=0.25) +
  geom_path(data=galellip, aes(x=x, y=y), color='black', linetype=1, size=0.25) +
  geom_point(aes(x=galcenx, y=galceny), color='black', shape=5) + 
  scale_x_continuous(limits=c(minx,maxx)) + 
  scale_y_continuous(limits=c(miny,maxy))

g4 <- g4 + theme_bw() +
  xlab("x") +
  ylab("y") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=5), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=5)) +
  annotate("text", x=Inf,y=-Inf
           , label =paste("B-R <=",max(db$BR)),size = 2, hjust = 1.1, vjust = -1.0)
  

g5 <- g5 + theme_bw() +
  xlab("x") +
  ylab("y") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=5), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=5)) +
  annotate("text", x=Inf,y=-Inf
           , label =paste(max(db$BR),"< B-R <=",max(dg$BR)),size = 2, hjust = 1.1, vjust = -1.0)

g6 <- g6 + theme_bw() +
  xlab("x") +
  ylab("y") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=5), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=5)) +
  annotate("text", x=Inf,y=-Inf
           , label =paste("B-R >",max(dg$BR)),size = 2, hjust = 1.1, vjust = -1.0)


lay <- rbind(c(1,1,1,2,2,2),c(3,3,3,4,4,4),c(5,5,6,6,7,7),c(5,5,6,6,7,7))


gplot <- arrangeGrob(g1, g3a, g2, g3b, g4, g5, g6, layout_matrix= lay, respect=TRUE)
ggsave(outfile, gplot, width=8, height=5.5)
