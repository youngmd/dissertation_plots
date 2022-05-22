library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(gtable)
library(mixtools)
library(nortest)
library(ADGofTest)
library(fitdistrplus)
# Read in the data

options(warn=-1)

lratioP <- function(crit, idf) {
  x <- crit
  an1 <- idf
  f <- x/an1
  an2 <- 1e10
  ff <- f
  prob <- 1.0

  if (f < 1.0) {
    ff <- 1.0 / f
    temp <- an1
    an1 <- an2
    an2 <- temp
  }

  a1 <- 2.0 / an1 / 9.0
  a2 <- 2.0 / an2 /9.0
  z <- abs(((1.0 - a2) * ff^0.3333 - 1.0 + a1)/sqrt(a2 * ff^0.666666 + a1))

  if (an2 <= 3.0) { z <- z*(1.0 + 0.08 * z^4 / an2^3)}

  fz <- exp(-z*z/2.0)*0.3989423
  w <-1.0 / (1.0 + z * 0.2316419)
  prob <- fz * w * ((((1.332074*w-1.821256)*w+1.781478)*w-0.3565638)*w + 0.3193815)

  if (f < 1.0) { prob <- 1.0 - prob}
  prob
}

mixture_crossing <- function(foo){
  N <- length(foo$mu) - 1
  crossings <- rep(0,N)
  for (i in 1:N) {
    crossing <- uniroot(function(x) dnorm(x, mean=foo$mu[i], sd=foo$sigma[i]) - dnorm(x, mean=foo$mu[i+1], sd=foo$sigma[i+1]), c(foo$mu[i],foo$mu[i+1]))
    crossings[i] <- crossing$root
  }
  return(crossings)
}

parseGccs <- function(datapath, galaxy){
  cnames1 <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr","raddist")

  cnames2 <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr")

  gccdata <- tryCatch({
      read.table(datapath, col.names=cnames1)
    },
    warning = function(w){
      print()
    },
    error = function(err){
      return(read.table(datapath, col.names=cnames2))
    }
  )
  gccdata$bmr <- gccdata$bmv + gccdata$vmr
  gccdata$galname <- galaxy
  return(data.frame(bmr = gccdata$bmr, raddist = gccdata$raddist, galname= gccdata$galname))
}

mixture_crossing <- function(foo){
  N <- length(foo$mu) - 1
  crossings <- rep(0,N)
  for (i in 1:N) {
    crossing <- uniroot(function(x) dnorm(x, mean=foo$mu[i], sd=foo$sigma[i]) - dnorm(x, mean=foo$mu[i+1], sd=foo$sigma[i+1]), c(foo$mu[i],foo$mu[i+1]))
    crossings[i] <- crossing$root
  }
  return(crossings)
}

fitmm <- function(gccdata, galname, binw){

  mu2 = c(1.2, 1.4)
  l2 = c(0.5, 0.5)
  sig2 = c(0.1, 0.1)

  m1 <- fitdist(gccdata$bmr, "norm")
  m2 <- normalmixEM(gccdata$bmr, mu=mu2, maxit=10000, lambda=l2, sigma=sig2, k=2, mean.constr=c(NA, NA), arbmean=TRUE, arbvar=TRUE)

  ll1 <- logLik(glm(gccdata$bmr~1))[1]

  m2$crossing <- mixture_crossing(m2)
  foo <- m2

  m2f_lr <- -2.0 * (ll1 - m2$loglik)
  m2f_p <- lratioP(m2f_lr, 4)
  
  print(nrow(gccdata))
  fblue <- subset(gccdata, bmr < 1.23)
  print(nrow(fblue) / nrow(gccdata))
  # binw <- 0.05
  print(galname)
  summary(m2)
  print(m2$crossing)
  dist <- data.frame( BmR = seq(from = 0.7, to = 1.9, length.out = 500))

  dist$B <- dnorm(dist$BmR,mean=m2$mu[1], sd=m2$sigma[1]) * (length(gccdata$bmr) * m2$lambda[1] * binw)
  dist$R <- dnorm(dist$BmR,mean=m2$mu[2], sd=m2$sigma[2]) * (length(gccdata$bmr) * m2$lambda[2] * binw)
  dist$galname <- galname
  return(dist)
}

galname <- "NGC 4382"
datapath <- "~/gcc/n4382/gcfinder3/gc_cand_calib.90.out"
datapath2 <- "~/gcc/n4382/gcfinder3/keep.out"

gcc90data <- parseGccs(datapath, galname)
gcc90data <- subset(gcc90data, bmr < 1.9)
x <- gccdata$bmr
binw <- 3 * IQR(x) / length(x)^(1/3)

fit90dist <- fitmm(gcc90data, galname, binw)
gccdata <- parseGccs(datapath2, galname)

galname <- "NGC 5846"
datapath <- "~/gcc/n5846/gcfinder3/gc_cand_calib.90.out"
datapath2 <- "~/gcc/n5846/gcfinder3/keep.out"

tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data <- rbind(gcc90data, tmp90)
fit90dist <- rbind(fit90dist, fitmm(tmp90, galname, binw))
gccdata <- rbind(gccdata, parseGccs(datapath2, galname))


galname <- "NGC 4649"
datapath <- "~/gcc/n4649/gcfinder3.n4649/gc_cand_calib.90.out"
datapath2 <- "~/gcc/n4649/gcfinder3.n4649/keep.out"

tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data <- rbind(gcc90data, tmp90)
fit90dist <- rbind(fit90dist, fitmm(tmp90, galname, binw))
gccdata <- rbind(gccdata, parseGccs(datapath2, galname))


galname <- "NGC 4621"
datapath <- "~/gcc/n4649/gcfinder3.n4621/gc_cand_calib.90.out"
datapath2 <- "~/gcc/n4649/gcfinder3.n4621/keep.out"

tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data <- rbind(gcc90data, tmp90)
fit90dist <- rbind(fit90dist, fitmm(tmp90, galname, binw))
gccdata <- rbind(gccdata, parseGccs(datapath2, galname))



gccdata <- subset(gccdata, bmr < 1.9)
gccdata <- subset(gccdata, bmr > 0.8)
gcc90data <- subset(gcc90data, bmr < 1.9)
gcc90data <- subset(gcc90data, bmr > 0.8)

gccdata$gal_s = factor(gccdata$galname, levels=c('NGC 4649','NGC 5846','NGC 4382', 'NGC 4621'))
gcc90data$gal_s = factor(gcc90data$galname, levels=c('NGC 4649','NGC 5846','NGC 4382', 'NGC 4621'))
fit90dist$gal_s = factor(fit90dist$galname, levels=c('NGC 4649','NGC 5846','NGC 4382', 'NGC 4621'))

g1 <- ggplot(data=fit90dist) + 
    geom_histogram(data=gcc90data, aes(x=bmr), fill="grey", binwidth=binw, alpha=0.8) +
    geom_histogram(data=gccdata, aes(x=bmr), binwidth=binw, alpha=0.3) +
    geom_path(aes( BmR, B), linetype=2, color='blue', size=0.3) +
    geom_path(aes( BmR, R), linetype=3, color='red', size= 0.3) +
    geom_blank(data=gccdata, aes(x=bmr, y=180), stat="identity") +
    facet_wrap(~gal_s,nrow=2) +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab(expression('Number')) + 
    xlab(expression('B - R')) + 
    scale_y_continuous(expand= c(0,0)) 


ggsave("gmm_mosaic.pdf", g1, width=8, height=8)


galname <- "NGC 1023"
datapath <- "~/gcc/n1023/gc_cand_calib.90.out"
datapath2 <- "~/gcc/n1023/keep.out"

x <- gccdata2$bmr
binw <- 3 * IQR(x) / length(x)^(1/3)

tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data2 <- tmp90
fit90dist2 <- fitmm(tmp90, galname, binw)
gccdata2 <- parseGccs(datapath2, galname)


galname <- "NGC 7332"
datapath <- "~/gcc/n7332/gc_cand_calib.90.dat"
datapath2 <- "~/gcc/n7332/keep.out"

tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data2 <- rbind(gcc90data2, tmp90)
fit90dist2 <- rbind(fit90dist2, fitmm(tmp90, galname, binw))
gccdata2 <- rbind(gccdata2, parseGccs(datapath2, galname))


galname <- "NGC 4013"
datapath <- "~/gcc/n4013/gc_cand_calib.90.dat"
datapath2 <- "~/gcc/n4013/keep.out"


tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data2 <- rbind(gcc90data2, tmp90)
fit90dist2 <- rbind(fit90dist2, fitmm(tmp90, galname, binw))
gccdata2 <- rbind(gccdata2, parseGccs(datapath2, galname))

galname <- "NGC 7339"
datapath <- "~/gcc/n7339/gc_cand_calib.90.dat"
datapath2 <- "~/gcc/n7339/keep.out"


tmp90 <- parseGccs(datapath, galname)
tmp90 <- subset(tmp90, bmr < 1.9)
gcc90data2 <- rbind(gcc90data2, tmp90)
fit90dist2 <- rbind(fit90dist2, fitmm(tmp90, galname, binw))
gccdata2 <- rbind(gccdata2, parseGccs(datapath2, galname))



gccdata2$gal_s = factor(gccdata2$galname, levels=c('NGC 1023', 'NGC 4013', 'NGC 7332', 'NGC 7339'))
gcc90data2$gal_s = factor(gcc90data2$galname, levels=c('NGC 1023', 'NGC 4013', 'NGC 7332', 'NGC 7339'))
fit90dist2$gal_s = factor(fit90dist2$galname, levels=c('NGC 1023', 'NGC 4013', 'NGC 7332', 'NGC 7339'))

gccdata2 <- subset(gccdata2, bmr < 1.9)
gccdata2 <- subset(gccdata2, bmr > 0.8)

g2 <- ggplot(fit90dist2) + 
    geom_histogram(data=gcc90data2, aes(x=bmr), fill="grey", binwidth=binw, alpha=0.8) +
    geom_histogram(data=gccdata2, aes(x=bmr), binwidth=binw, alpha=0.3) +
    geom_path(aes( BmR, B), linetype=2, color='blue', size=0.3) +
    geom_path(aes( BmR, R), linetype=3, color='red', size= 0.3) +
    geom_blank(data=gccdata2, aes(x=bmr, y=60), stat="identity") +
    facet_wrap(~gal_s,nrow=2) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab(expression('Number')) + 
    xlab(expression('B - R')) + 
    scale_y_continuous(expand= c(0,0)) 


ggsave("gmm_wiyn.pdf", g2, width=8, height=8)

#plot(mirror.ticks(g))
# gplot <- arrangeGrob(g1, nrow=3, ncol=1)
# ggsave(outfile, gplot, width=8, height=8)
#ggsave(outfile, g, width=7, height=7)
#embed_fonts('tmp.pdf', outfile=outfile)

# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#     r = diameter / 2
#     tt <- seq(0,2*pi,length.out = npoints)
#     xx <- center[1] + r * cos(tt)
#     yy <- center[2] + r * sin(tt)
#     return(data.frame(x = xx, y = yy))
# }

# dat <- circleFun(c(xcen,ycen),extent,npoints = 100)

# bgccs <- gccdata[gccdata$bmr < m2join$crossing,]
# rgccs <- gccdata[gccdata$bmr >= m2join$crossing,]

# s1 <- ggplot(data=bgccs) +
#       geom_point(aes(x=bx, y=by), color="blue") +
#       geom_path(data=dat, aes(x,y), linetype=3) +
#       theme_bw() + 
#       theme(legend.position="none") + 
#       theme(text=element_text(family="CM Roman")) +
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(), 
#             panel.border = element_blank(), 
#             axis.line = element_line(colour = "black"),
#             axis.ticks.length=unit(-0.25, "cm"),
#             axis.text.y=element_text(margin=margin(0,15,0,0)), 
#             axis.text.x=element_text(margin=margin(15,0,0,0))) +
#       # annotate("text", x = Inf, y = Inf, label = g4label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
#       ylab(expression('Y position (pixels)')) + 
#       xlab(expression('X position (pixels)'))

# s2 <- ggplot(data=rgccs) +
#       geom_point(aes(x=bx, y=by), color="red") +
#       geom_path(data=dat, aes(x,y), linetype=3) +
#       theme_bw() + 
#       theme(legend.position="none") + 
#       theme(text=element_text(family="CM Roman")) +
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(), 
#             panel.border = element_blank(), 
#             axis.line = element_line(colour = "black"),
#             axis.ticks.length=unit(-0.25, "cm"),
#             axis.text.y=element_text(margin=margin(0,15,0,0)), 
#             axis.text.x=element_text(margin=margin(15,0,0,0))) +
#       # annotate("text", x = Inf, y = Inf, label = g4label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
#       ylab(expression('Y position (pixels)')) + 
#       xlab(expression('X position (pixels)'))

# splot <- arrangeGrob(s1, s2, s2, s1, nrow=2, ncol=2)
# ggsave(outfile2, splot, width=8, height=8)
# x <- bmr

# N <- length(x)
# ks.boot <- rep(0,1000)
# for (i in 1:1000) {
#   z <- rbinom(N, 1, foo$lambda[1])
#   x.b <- z*rnorm(N, foo$mu[1], foo$sigma[1]) + (1-z)*rnorm(N, foo$mu[2], foo$sigma[2])
#   foo.b <- normalmixEM(x.b, maxit=10000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma)
#   ks.boot[i] <- ks.test(x.b, pmnorm, mu=foo.b$mu, sigma=foo.b$sigma, pmix=foo.b$lambda)$statistic
# }

# mean(test$statistic <= ks.boot)