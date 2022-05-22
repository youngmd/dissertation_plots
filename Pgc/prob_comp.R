library(ggplot2)
data <- read.table('n4594_hrv.tsv', header=TRUE)
data$vi <- abs(data$HRV - 1091)
data$mn <- data$vi^2 * data$Rad
n <- nrow(data) - 219
datasort <- data[order(data$mn),]
gccs <- datasort[1:n,]
cutoff <- max(gccs$mn)
data$isGCC[data$mn <= cutoff] <- TRUE
data$isGCC[data$mn > cutoff] <- FALSE


contams <- read.csv('n4594_contam.csv', header=TRUE)
z <- contams$V
binw <- 4 * IQR(z) / length(z)^(1/3)

nbins <- round((max(z) - min(z)) / binw)
contamdist <- hist(z, breaks='scott')
# ggplot(data=data, aes(x=Rad, y=HRV))+geom_point()+geom_point(data=gccs, aes(x=Rad, y=HRV), color='blue')+scale_y_continuous(limits=c(-1000,3000))

gclfdat <- read.table('n4594_gclf.dat', header=TRUE)

dev_a0 <- 3.76
da0_err <- 0.11
da1_err <- 0.08
dev_a1 <- -2.11
pow_a0 <- 1.87
pa0_err <- 0.05
pow_a1 <- -1.85
pa1_err <- 0.07
contam <- 0.54 # not done asymptotically, stars+galaxies = 0.26 + 0.28
contam_err <- 0.10
gclf <- 0.83
gclf_err <- 0.02
extent <- 19.0

gclf_peak <- 22.6

x1 <- rnorm(1e5, gclf_peak, 1.2)
x2 <- rnorm(1e5, gclf_peak, 1.3)
x3 <- rnorm(1e5, gclf_peak, 1.4)

x <- c(x1, x2, x3)

devPAlt <- function(r, v){
  da0 <- dev_a0
  da1 <- dev_a1
  c <- contam
  lsigma = da0 + da1*(r^0.25)
  sigma = 10^lsigma
  probs = vector(,length(r))
  for(j in 1:length(sigma)){
    for(i in 1:length(contamdist$breaks)){
       if(contamdist$breaks[i+1] > v[j]){
         fraccontam = contamdist$density[i]
         lwr = contamdist$breaks[i]
         upr = contamdist$breaks[i+1]
         nx <- x[x < upr & x > lwr]
         compl <- approx(gclfdat$Vmag, gclfdat$compl, v[j])
         print(v[j])
         print(compl$y)
         fracgcc <- ( length(nx)/ length(x) )
         siggcc <- sigma[j] * fracgcc * compl$y
         sigcontam <- c * fraccontam
         prob <- siggcc / (siggcc + sigcontam)
         probs[j] <- prob
         break
       }
    }
  }
  return(probs)
}


powP <- function(r,err=FALSE){
    pa0 <- pow_a0
    pa1 <- pow_a1
    g <- gclf
    c <- contam
    if(err == 'up'){
        pa0 <- pa0 + pa0_err
        pa1 <- pa1 + pa1_err
        g <- g + gclf_err
        c <- c - contam_err
    }
    if(err == 'down'){
        pa0 <- pa0 - pa0_err
        pa1 <- pa1 - pa1_err
        g <- g - gclf_err
        c <- c + contam_err
    }
    lsigma = pa0 + pa1*log10(r)
    sigma = 10^lsigma
    corr_sigma = sigma * g
    prob = corr_sigma / (corr_sigma + c)
    return(prob)
}

devP <- function(r,err=FALSE){
    da0 <- dev_a0
    da1 <- dev_a1
    g <- gclf
    c <- contam
    if(err == 'up'){
        da0 <- da0 + da0_err
        da1 <- da1 + da1_err
        g <- g + gclf_err
        c <- c - contam_err
    }
    if(err == 'down'){
        da0 <- da0 - da0_err
        da1 <- da1 - da1_err
        g <- g - gclf_err
        c <- c + contam_err
    }
    lsigma = da0 + da1*(r^0.25)
    sigma = 10^lsigma
    corr_sigma = sigma * g
    prob = corr_sigma / (corr_sigma + c)
    return(prob)            
}


data$gccPdev <- devPAlt(data$Rad, data$Vmag)
data$gccPdevold <- devP(data$Rad)
# data$gccPdev_up <- devP(data$Rad, 'up')
# data$gccPdev_down <- devP(data$Rad, 'down')
# 
# data$gccPpow <- powP(data$Rad)
# data$gccPpow_up <- powP(data$Rad, 'up')
# data$gccPpow_down <- powP(data$Rad, 'down')
datasort <- data[order(data$Rad),]

rz <- subset(datasort, r_Name =='RZ')
rz$gccPdevco <- cumsum(rz$gccPdevold)
# rz$gccPpowc <- cumsum(rz$gccPpow)
# rz$gccPpowdc <- cumsum(rz$gccPpow_down)
# rz$gccPpowuc <- cumsum(rz$gccPpow_up)

rz$gccPdevc <- cumsum(rz$gccPdev)
#rz$gccPdevdc <- cumsum(rz$gccPdev_down)
#rz$gccPdevuc <- cumsum(rz$gccPdev_up)

rz$isGCCc <- cumsum(rz$isGCC)

#cols <- c("Cumulative Prediction"="#f04546","Cumulative Observed"="#3591d1","Prediction Interval"="grey90")
cols <- c("Cumulative Prediction"="#f04546","Cumulative Observed"="#3591d1","Cumulative Prediction (old)"="#25A586")

cols2 <- c("GC"="#3591d1","Not GC"="#f04546")

g1 <- ggplot(data=rz, aes(x=Rad, y=isGCCc)) + 
    #geom_ribbon(aes(ymin=gccPdevdc, ymax=gccPdevuc, fill='Prediction Interval')) + 
    geom_line(aes(y=gccPdevc, color="Cumulative Prediction"), linetype=1, size=0.3) +
    #geom_line(aes(y=gccPdevco, color="Cumulative Prediction (old)"), linetype=1, size=0.3) + 
    geom_line(aes(y=isGCCc, color="Cumulative Observed"), linetype=1, size=0.3) +
    scale_y_continuous(limits=c(0,500)) + 
    ylab("Total") +
    xlab("Projected Radial Distance (arcmin)") +
    ggtitle(expression(Cumulative~P[GC]~compared~to~spectroscopic~results~(NGC4594))) +
    scale_colour_manual(name="", values=cols) + 
    scale_linetype_manual(values=c(2,3)) + 
    scale_fill_manual(name="", values=cols) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          legend.position=c(0.2,0.7),
          legend.key = element_blank(),
 #         panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank())

g2 <- ggplot(data=rz, aes(x=Rad, y=HRV)) + 
        geom_point(aes(color=isGCC, shape=isGCC), size=1) +
        geom_hline(aes(yintercept=1099.5), linetype=2, size=0.2, color='red4') +
        scale_y_continuous(limits=c(-1000,3000)) + 
        ylab(expression(v[helio]~(km~s^-1))) +
        xlab("Projected Radial Distance (arcmin)") +
        theme_bw() + 
        theme(legend.position="none") + 
        theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          legend.key = element_blank(),
#          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0)))

g3 <- rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last")
# g3 <- grid.arrange(g1,g2, ncol=1, heights=c(2,1))
ggsave("Pgc2_alt_4594.pdf", g3, width=8, height=8)