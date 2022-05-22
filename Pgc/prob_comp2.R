library(ggplot2)
data <- read.table('n3379.dat', header=TRUE)
# data$vi <- abs(data$HRV - 1091)
# data$mn <- data$vi^2 * data$Rad
# n <- nrow(data) - 219
# datasort <- data[order(data$mn),]
# gccs <- datasort[1:n,]
# cutoff <- max(gccs$mn)
data$isGCC[data$ngc == 'Y'] <- TRUE
data$isGCC[data$ngc == 'N'] <- FALSE
data$isGCC[data$v > 578] <- TRUE
# ggplot(data=data, aes(x=Rad, y=HRV))+geom_point()+geom_point(data=gccs, aes(x=Rad, y=HRV), color='blue')+scale_y_continuous(limits=c(-1000,3000))

dev_a0 <- 2.33
da0_err <- 0.22
da1_err <- 0.15
dev_a1 <- -1.62
pow_a0 <- 0.83
pa0_err <- 0.09
pow_a1 <- -1.41
pa1_err <- 0.13
contam <- 0.29 # not done asymptotically, stars+galaxies = 0.08 + 0.11
contam_err <- 0.02
gclf <- 0.46
gclf_err <- 0.01
extent <- 11.0

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


data$gccPdev <- devP(data$Rad)
data$gccPdev_up <- devP(data$Rad, 'up')
data$gccPdev_down <- devP(data$Rad, 'down')

data$gccPpow <- powP(data$Rad)
data$gccPpow_up <- powP(data$Rad, 'up')
data$gccPpow_down <- powP(data$Rad, 'down')
datasort <- data[order(data$Rad),]

rz <- datasort
rz$gccPpowc <- cumsum(rz$gccPpow)
rz$gccPpowdc <- cumsum(rz$gccPpow_down)
rz$gccPpowuc <- cumsum(rz$gccPpow_up)

rz$gccPdevc <- cumsum(rz$gccPdev)
rz$gccPdevdc <- cumsum(rz$gccPdev_down)
rz$gccPdevuc <- cumsum(rz$gccPdev_up)

rz$isGCCc <- cumsum(rz$isGCC)

cols <- c("Cumulative Prediction"="#f04546","Cumulative Observed"="#3591d1","Prediction Interval"="grey90")

cols2 <- c("GC"="#3591d1","Not GC"="#f04546")

g1 <- ggplot(data=rz, aes(x=Rad, y=isGCCc)) + 
    geom_ribbon(aes(ymin=gccPdevdc, ymax=gccPdevuc, fill='Prediction Interval')) + 
    geom_line(aes(y=gccPdevc, color="Cumulative Prediction"), linetype=1, size=0.3) + 
    geom_line(aes(y=isGCCc, color="Cumulative Observed"), linetype=1, size=0.3) +
    scale_y_continuous(limits=c(0,150)) + 
    ylab("Total") +
    xlab("Projected Radial Distance (arcmin)") +
    ggtitle(expression(Cumulative~P[GC]~compared~to~spectroscopic~results~(NGC3379))) +
    scale_colour_manual(name="", values=cols) + 
    scale_linetype_manual(values=c(2,3)) + 
    scale_fill_manual(name="", values=cols) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          legend.position=c(0.2,0.7),
          legend.key = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank())

g2 <- ggplot(data=rz, aes(x=Rad, y=v)) + 
        geom_point(aes(color=isGCC, shape=isGCC), size=1) +
        geom_hline(aes(yintercept=916.7), linetype=2, size=0.2, color='red4') +
        scale_y_continuous(limits=c(-300,1200)) + 
        ylab(expression(v[helio]~(km~s^-1))) +
        xlab("Projected Radial Distance (arcmin)") +
        theme_bw() + 
        theme(legend.position="none") + 
        theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          legend.key = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0)))

g3 <- rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last")
# g3 <- grid.arrange(g1,g2, ncol=1, heights=c(2,1))
ggsave("Pgc3379.pdf", g3, width=8, height=8)