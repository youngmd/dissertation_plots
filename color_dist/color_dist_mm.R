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
#  PROB (CRIT,IDF)
# C     For accuracy see Golden, Weiss and Davis (1968)
# C     Educ. Psycol. Measurement Vol 28, 163-165.
# C     Reference R.G. Davies. Computer Programming in Quantitative
# C     Biology.
# C
#       X = CRIT
#       AN1 = IDF
#       F=X/AN1
#       AN2=1.0E10
#       FF=F
#       PROB=1.0
#       IF (AN1*AN2*F.EQ.0.0) RETURN
# C
# C     Take reciprocal if F les than 1.0
# C
#       IF (F.GE.1.0) GO TO 6
#       FF=1.0/F
#       TEMP=AN1
#       AN1=AN2
#       AN2=TEMP
# C
# C     Normalise Variance Ratio
# C
# 6     A1=2.0/AN1/9.0
#       A2=2.0/AN2/9.0
#       Z=ABS(((1.0-A2)*FF**0.3333333-1.0+A1)/SQRT(A2*FF**0.6666666+A1))
#       IF (AN2.LE.3.0) Z=Z*(1.0+0.08*Z**4/AN2**3)
# C
# C     Compute probability
# C
# 7     FZ = EXP(-Z*Z/2.0)*0.3989423
#       W=1.0/(1.0+Z*0.2316419)
#       PROB = FZ*W*((((1.332074*W-1.821256)*W+1.781478)*W-0.3565638)*W
#      1+0.3193815)
#       IF (F.LT.1.0) PROB=1.0-PROB
#       RETURN
#       END

#2 mix
pmnorm2 <- function(x, mu, sigma, pmix) {
  pmix[1]*pnorm(x,mu[1],sigma[1]) + (1-pmix[1])*pnorm(x,mu[2],sigma[2])
}

#3 mix
pmnorm3 <- function(x, mu, sigma, pmix) {
  pmix[1]*pnorm(x,mu[1],sigma[1]) + (pmix[2])*pnorm(x,mu[2],sigma[2]) + (pmix[3])*pnorm(x,mu[3],sigma[3])
}

pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  library(extrafont)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# bootstrapKS2 <- function(x, foo){
#   N <- length(x)
#   ks.boot <- rep(0,100)
#   for (i in 1:100) {
#     z <- rbinom(N, 1, foo$lambda[1])
#     x.b <- z*rnorm(N, foo$mu[1], foo$sigma[1]) + (1-z)*rnorm(N, foo$mu[2], foo$sigma[2])
#     foo.b <- normalmixEM(x.b, maxit=10000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma)
#     ks.boot[i] <- ks.test(x.b, pmnorm2, mu=foo.b$mu, sigma=foo.b$sigma, pmix=foo.b$lambda)$statistic
#   }

#   mean(test$statistic <= ks.boot)
# }

ADiters <- 100

# bootstrapKS2 <- function(x, foo, arbvar=TRUE, burnin=25){
#   test <- ks.test(x, pnormmix, mixture=foo)
#   if(test$p.value < 0.05){  #if even the naive test rejects the null, no point doing the rest
#     return(test$p.value)
#   }
#   N <- length(x)
#   ks.boot <- rep(0,KSiters)
#   for (i in 1:KSiters) {
#     z <- rbinom(N, 1, foo$lambda[1])
#     x.b <- z*rnorm(N, foo$mu[1], foo$sigma[1]) + (1-z)*rnorm(N, foo$mu[2], foo$sigma[2])
#     foo.b <- normalmixEM(x.b, maxit=10000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma, arbvar=arbvar)
#     ks.boot[i] <- ks.test(x.b, pnormmix, mixture=foo.b)$p.value
#     if(i > burnin){
#       #print(mean(ks.boot))
#       if(sd(ks.boot) < 0.005){
#         cat("Converged in loop ", i)
#         break
#       }
#     }
#   }

#   mean(ks.boot)
# }

bootstrapAD <- function(x, foo){
  test <- ADGofTest::ad.test(x, pnorm, mean=foo$estimate["mean"], sd=foo$estimate["sd"])
  if(test$p.value < 0.05){  #if even the naive test rejects the null, no point doing the rest
    return(test$p.value)
  }
  N <- length(x)
  ad.boot <- rep(0,ADiters)
  for (i in 1:ADiters) {
    x.b <- rnorm(N, foo$estimate["mean"], foo$estimate["sd"])
    foo.b <- fitdist(x.b, "norm")
    ad.boot[i] <- ADGofTest::ad.test(x.b, pnorm, mean=foo.b$estimate["mean"], sd=foo.b$estimate["sd"])$statistic
  }

  mean(test$statistic <= ad.boot)
}

bootstrapAD2 <- function(x, foo, arbvar=TRUE){
  test <- ADGofTest::ad.test(x, pnormmix, mixture=foo)
  if(test$p.value < 0.05){  #if even the naive test rejects the null, no point doing the rest
    return(test$p.value)
  }
  N <- length(x)
  ad.boot <- rep(0,ADiters)
  for (i in 1:ADiters) {
    z <- rbinom(N, 1, foo$lambda[1])
    x.b <- z*rnorm(N, foo$mu[1], foo$sigma[1]) + (1-z)*rnorm(N, foo$mu[2], foo$sigma[2])
    foo.b <- normalmixEM(x.b, maxit=20000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma, arbvar=arbvar)
    ad.boot[i] <- ADGofTest::ad.test(x.b, pnormmix, mixture=foo.b)$statistic
  }

  mean(test$statistic <= ad.boot)
}


bootstrapAD3 <- function(x, foo, arbvar=TRUE){
  test <- ADGofTest::ad.test(x, pnormmix, mixture=foo)
  if(test$p.value < 0.05){ #if even the naive test rejects the null, no point doing the rest
    return(test$p.value)
  }
  N <- length(x)
  ad.boot <- rep(0,ADiters)
  for (i in 1:ADiters) {
    z <- rmultinom(N, 2, foo$lambda)
    z1 <- z
    z1[which(z==0)] = 1
    z1[which(z!=0)] = 0
    z2 <- z
    z2[which(z!=1)] = 0
    z3 <- z
    z3[which(z==2)] = 1
    z3[which(z!=2)] = 0
    x.b <- z1*rnorm(N, foo$mu[1], foo$sigma[1]) + z2*rnorm(N, foo$mu[2], foo$sigma[2]) + z3*rnorm(N, foo$mu[3], foo$sigma[3])
    foo.b <- normalmixEM(x.b, maxit=10000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma, arbvar=arbvar)
    ad.boot[i] <- ADGofTest::ad.test(x.b, pnormmix, mixture=foo.b)$statistic
  }

  mean(test$statistic <= ad.boot)
}
# bootstrapKS3 <- function(x, foo, arbvar=TRUE, burnin=25){
#   test <- ks.test(x, pnormmix, mixture=foo)
#   if(test$p.value < 0.05){ #if even the naive test rejects the null, no point doing the rest
#     return(test$p.value)
#   }
#   N <- length(x)
#   ks.boot <- rep(0,KSiters)
#   for (i in 1:KSiters) {
#     z <- rmultinom(N, 2, foo$lambda)
#     z1 <- z
#     z1[which(z==0)] = 1
#     z1[which(z!=0)] = 0
#     z2 <- z
#     z2[which(z!=1)] = 0
#     z3 <- z
#     z3[which(z==2)] = 1
#     z3[which(z!=2)] = 0
#     x.b <- z1*rnorm(N, foo$mu[1], foo$sigma[1]) + z2*rnorm(N, foo$mu[2], foo$sigma[2]) + z3*rnorm(N, foo$mu[3], foo$sigma[3])
#     foo.b <- normalmixEM(x.b, maxit=10000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma, arbvar=arbvar)
#     ks.boot[i] <- ks.test(x.b, pnormmix, mixture=foo.b)$p.value
#     if(i > burnin){
#       #print(mean(ks.boot))
#       if(sd(ks.boot) < 0.005){
#         cat("Converged in loop ", i)
#         break
#       }
#     }
#   }

#   mean(ks.boot)
# }

mixture_crossing <- function(foo){
  N <- length(foo$mu) - 1
  crossings <- rep(0,N)
  for (i in 1:N) {
    crossing <- uniroot(function(x) dnorm(x, mean=foo$mu[i], sd=foo$sigma[i]) - dnorm(x, mean=foo$mu[i+1], sd=foo$sigma[i+1]), c(foo$mu[i],foo$mu[i+1]))
    crossings[i] <- crossing$root
  }
  return(crossings)
}

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
gccdata <- subset(gccdata, bmr < 1.9)
gccdata <- subset(gccdata, bmr > 0.8)
bmr <- gccdata$bmv + gccdata$vmr

mu2 = c(1.1, 1.4)
mu3 = c(1.0, 1.3, 1.5)
l2 = c(0.5, 0.5)
l3 = c(0.3, 0.4, 0.3)
sig2 = c(0.1, 0.1)
sig3 = c(0.1, 0.1, 0.1)

m1 <- fitdist(bmr, "norm")
m2 <- normalmixEM(bmr, mu=mu2, maxit=10000, lambda=l2, sigma=sig2, k=2, mean.constr=c(NA, NA), arbmean=TRUE, arbvar=TRUE)
m3 <- normalmixEM(bmr, mu=mu3, maxit=10000, lambda=l3, sigma=sig3, k=3, arbmean=TRUE, arbvar=TRUE)

ll1 <- logLik(glm(bmr~1))[1]

m2$crossing <- mixture_crossing(m2)
m3$crossing <- mixture_crossing(m3)
foo <- m2

m2f_lr <- -2.0 * (ll1 - m2$loglik)
m2f_p <- lratioP(m2f_lr, 4)

m3f_lr <- -2.0 * (ll1 - m3$loglik)
m3f_p <- lratioP(m3f_lr, 8)

x <- gccdata$bmr
binw <- 2 * IQR(x) / length(x)^(1/3)
# binw <- 0.05
dist <- data.frame( BmR = seq(from = 0.7, to = 1.9, length.out = 500))

# dist$dist1.3 <- dnorm(dist$Vmag,mean=23.93, sd=1.3) * 565.871972
# dist$dist1.4 <- dnorm(dist$Vmag,mean=23.93, sd=1.4) * 541.143759

print("Testing unimodal case")
log1 <- capture.output(m1ad <- bootstrapAD(bmr, m1))

print("Bootstrapping A-D statistic")
ADiters

log2 <- capture.output(m2ad <- bootstrapAD2(bmr, m2))
log3 <- capture.output(m3ad <- bootstrapAD3(bmr, m3))

m1ad
m2ad
m3ad

dist$U <- dnorm(dist$BmR,mean=m1$estimate["mean"], sd=m1$estimate["sd"]) * (length(x) * binw)
dist$B <- dnorm(dist$BmR,mean=m2$mu[1], sd=m2$sigma[1]) * (length(bmr) * m2$lambda[1] * binw)
dist$R <- dnorm(dist$BmR,mean=m2$mu[2], sd=m2$sigma[2]) * (length(bmr) * m2$lambda[2] * binw)

g0label <- sprintf("Unimodal Gaussian\nll=%5.2f A-D=%4.3f", ll1, m1ad)
g1label <- sprintf("2 Gaussian Mixture\nllr=%5.2f p=%4.3f\nA-D=%4.3f\nMixing Fractions: %3.2f %3.2f\nCrossing Point: %4.3f ", m2f_lr, m2f_p, m2ad, m2$lambda[1], m2$lambda[2], m2$crossing)
g3label <- sprintf("3 Gaussian Mixture\nllr=%5.2f p=%4.3f\nA-D=%4.3f\nMixing Fractions: %3.2f %3.2f %3.2f\nCrossing Points: %4.3f %4.3f ", m3f_lr, m3f_p, m3ad, m3$lambda[1], m3$lambda[2], m3$lambda[3], m3$crossing[1],m3$crossing[2])


g0 <- ggplot() +
    geom_histogram(data=gccdata, aes(x=(bmv + vmr)), fill="grey", binwidth=binw) +
    # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
    geom_path(data = dist, aes( BmR, U), linetype=2, color='blue') +
    # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
    # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
    # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    annotate("text", x = Inf, y = Inf, label = g0label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
    ylab(expression('Number')) + 
    xlab(expression('B - R')) + 
    scale_y_continuous(expand= c(0,0)) 


g1 <- ggplot() + 
    geom_histogram(data=gccdata, aes(x=(bmv + vmr)), fill="grey", binwidth=binw) +
    # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
    geom_path(data = dist, aes( BmR, B), linetype=2, color='blue') +
    geom_path(data = dist, aes( BmR, R), linetype=3, color='red') +
    # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
    # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
    # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    annotate("text", x = Inf, y = Inf, label = g1label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
    ylab(expression('Number')) + 
    xlab(expression('B - R')) + 
    scale_y_continuous(expand= c(0,0)) 
    #annotate("text", x = 1.5, y = 0.2, label = galname,size = 6,family="CM Roman", hjust = 0, vjust = 0)


dist$B <- dnorm(dist$BmR,mean=m3$mu[1], sd=m3$sigma[1]) * (length(bmr) * m3$lambda[1] * binw)
dist$V <- dnorm(dist$BmR,mean=m3$mu[2], sd=m3$sigma[2]) * (length(bmr) * m3$lambda[2] * binw)
dist$R <- dnorm(dist$BmR,mean=m3$mu[3], sd=m3$sigma[3]) * (length(bmr) * m3$lambda[3] * binw)
crossing <- dist$BmR[which.min(abs(dist$distB - dist$distR))]

g3 <- ggplot() + 
    geom_histogram(data=gccdata, aes(x=(bmv + vmr)), fill="grey", binwidth=binw) +
    # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
    geom_path(data = dist, aes( BmR, B), linetype=2, color='blue') +
    geom_path(data = dist, aes( BmR, V), linetype=12, color='green4') +
    geom_path(data = dist, aes( BmR, R), linetype=3, color='red') +
    # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
    # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
    # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    annotate("text", x = Inf, y = Inf, label = g3label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
    ylab(expression('Number')) + 
    xlab(expression('B - R')) + 
    scale_y_continuous(expand= c(0,0)) 
    #annotate("text", x = 1.5, y = 0.2, label = galname,size = 6,family="CM Roman", hjust = 0, vjust = 0)

#plot(mirror.ticks(g))
gplot <- arrangeGrob(g0, g1, g3, nrow=3, ncol=1)
ggsave(outfile, gplot, width=8, height=8)
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