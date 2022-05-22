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
library(mclust)
# Read in the data

options(warn=-1)

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

ADiters <- 100

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
  ad.boot
  #mean(test$statistic <= ad.boot)
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
  ad.boot
  #mean(test$statistic <= ad.boot)
}

gcc90<- read.csv('gcc90.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

results <- read.table('mm_comp.dat', header=TRUE)
for(i in levels(gcc90$galaxy_name)){
  r <- subset(results, name == i)
  if(nrow(r) < 1){ next }
  gccdata <- subset(gcc90, galaxy_name == i)
  mm <- r[1,]$nmm
  if(mm == 1){ next }
  m <- Mclust(gccdata$BR, modelName = 'E', G=mm)

  if( mm == 2){
    m2 <- normalmixEM(gccdata$BR, mu=m$parameters$mean, maxit=10000, lambda=m$parameters$pro, sigma=sqrt(m$parameters$variance$sigmasq), k=2, mean.constr=c(NA, NA), arbmean=TRUE, arbvar=FALSE)
    log <- capture.output(ad <- bootstrapAD2(gccdata$BR, m2, arbvar=FALSE))
    ad.score <- ADGofTest::ad.test(gccdata$BR, pnormmix, m2)
  }

  if( mm == 3){
    m3 <- normalmixEM(gccdata$BR, mu=m$parameters$mean, maxit=10000, lambda=m$parameters$pro, sigma=sqrt(m$parameters$variance$sigmasq), k=3, mean.constr=c(NA, NA, NA), arbmean=TRUE, arbvar=FALSE)
    log <- capture.output(ad <- bootstrapAD3(gccdata$BR, m3, arbvar=FALSE))
    ad.score <- ADGofTest::ad.test(gccdata$BR, pnormmix, m3)
  }

  adp <- mean(ad.score$statistic <= ad)
  print(i)
  print(adp)
}


# mu2 = c(1.1, 1.4)
# mu3 = c(1.0, 1.3, 1.5)
# l2 = c(0.5, 0.5)
# l3 = c(0.3, 0.4, 0.3)
# sig2 = c(0.1, 0.1)
# sig3 = c(0.1, 0.1, 0.1)

# m1 <- fitdist(bmr, "norm")
# m2 <- normalmixEM(bmr, mu=mu2, maxit=10000, lambda=l2, sigma=sig2, k=2, mean.constr=c(NA, NA), arbmean=TRUE, arbvar=TRUE)
# m3 <- normalmixEM(bmr, mu=mu3, maxit=10000, lambda=l3, sigma=sig3, k=3, arbmean=TRUE, arbvar=TRUE)

# ll1 <- logLik(glm(bmr~1))[1]

# # m2$crossing <- mixture_crossing(m2)
# # m3$crossing <- mixture_crossing(m3)
# foo <- m2

# m2f_lr <- -2.0 * (ll1 - m2$loglik)
# m2f_p <- lratioP(m2f_lr, 4)

# m3f_lr <- -2.0 * (ll1 - m3$loglik)
# m3f_p <- lratioP(m3f_lr, 8)

# x <- gccdata$bmr
# binw <- 2 * IQR(x) / length(x)^(1/3)
# # binw <- 0.05
# dist <- data.frame( BmR = seq(from = 0.7, to = 1.9, length.out = 500))

# # dist$dist1.3 <- dnorm(dist$Vmag,mean=23.93, sd=1.3) * 565.871972
# # dist$dist1.4 <- dnorm(dist$Vmag,mean=23.93, sd=1.4) * 541.143759

# print("Testing unimodal case")
# log1 <- capture.output(m1ad <- bootstrapAD(bmr, m1))

# print("Bootstrapping A-D statistic")
# ADiters

# log2 <- capture.output(m2ad <- bootstrapAD2(bmr, m2))
# log3 <- capture.output(m3ad <- bootstrapAD3(bmr, m3))

# m1ad
# m2ad
# m3ad

# dist$U <- dnorm(dist$BmR,mean=m1$estimate["mean"], sd=m1$estimate["sd"]) * (length(x) * binw)
# dist$B <- dnorm(dist$BmR,mean=m2$mu[1], sd=m2$sigma[1]) * (length(bmr) * m2$lambda[1] * binw)
# dist$R <- dnorm(dist$BmR,mean=m2$mu[2], sd=m2$sigma[2]) * (length(bmr) * m2$lambda[2] * binw)

# g0label <- sprintf("Unimodal Gaussian\nll=%5.2f A-D=%4.3f", ll1, m1ad)
# g1label <- sprintf("2 Gaussian Mixture\nllr=%5.2f p=%4.3f\nA-D=%4.3f\nMixing Fractions: %3.2f %3.2f", m2f_lr, m2f_p, m2ad, m2$lambda[1], m2$lambda[2])
# g3label <- sprintf("3 Gaussian Mixture\nllr=%5.2f p=%4.3f\nA-D=%4.3f\nMixing Fractions: %3.2f %3.2f %3.2f", m3f_lr, m3f_p, m3ad, m3$lambda[1], m3$lambda[2], m3$lambda[3])


# g0 <- ggplot() +
#     geom_histogram(data=gccdata, aes(x=(bmr)), fill="grey", binwidth=binw) +
#     # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
#     geom_path(data = dist, aes( BmR, U), linetype=2, color='blue') +
#     # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
#     # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
#     # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
#     theme_bw() + 
#     theme(legend.position="none") + 
#     theme(text=element_text(family="CM Roman")) +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(), 
#           axis.line = element_line(colour = "black"),
#           axis.ticks.length=unit(-0.25, "cm"),
#           axis.text.y=element_text(margin=margin(0,15,0,0)), 
#           axis.text.x=element_text(margin=margin(15,0,0,0))) +
#     annotate("text", x = Inf, y = Inf, label = g0label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
#     ylab(expression('Number')) + 
#     xlab(expression('B - R')) + 
#     scale_y_continuous(expand= c(0,0)) 


# g1 <- ggplot() + 
#     geom_histogram(data=gccdata, aes(x=(bmr)), fill="grey", binwidth=binw) +
#     # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
#     geom_path(data = dist, aes( BmR, B), linetype=2, color='blue') +
#     geom_path(data = dist, aes( BmR, R), linetype=3, color='red') +
#     # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
#     # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
#     # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
#     theme_bw() + 
#     theme(legend.position="none") + 
#     theme(text=element_text(family="CM Roman")) +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(), 
#           axis.line = element_line(colour = "black"),
#           axis.ticks.length=unit(-0.25, "cm"),
#           axis.text.y=element_text(margin=margin(0,15,0,0)), 
#           axis.text.x=element_text(margin=margin(15,0,0,0))) +
#     annotate("text", x = Inf, y = Inf, label = g1label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
#     ylab(expression('Number')) + 
#     xlab(expression('B - R')) + 
#     scale_y_continuous(expand= c(0,0)) 
#     #annotate("text", x = 1.5, y = 0.2, label = galname,size = 6,family="CM Roman", hjust = 0, vjust = 0)


# dist$B <- dnorm(dist$BmR,mean=m3$mu[1], sd=m3$sigma[1]) * (length(bmr) * m3$lambda[1] * binw)
# dist$V <- dnorm(dist$BmR,mean=m3$mu[2], sd=m3$sigma[2]) * (length(bmr) * m3$lambda[2] * binw)
# dist$R <- dnorm(dist$BmR,mean=m3$mu[3], sd=m3$sigma[3]) * (length(bmr) * m3$lambda[3] * binw)
# crossing <- dist$BmR[which.min(abs(dist$distB - dist$distR))]

# g3 <- ggplot() + 
#     geom_histogram(data=gccdata, aes(x=(bmr)), fill="grey", binwidth=binw) +
#     # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
#     geom_path(data = dist, aes( BmR, B), linetype=2, color='blue') +
#     geom_path(data = dist, aes( BmR, V), linetype=12, color='green4') +
#     geom_path(data = dist, aes( BmR, R), linetype=3, color='red') +
#     # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
#     # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
#     # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
#     theme_bw() + 
#     theme(legend.position="none") + 
#     theme(text=element_text(family="CM Roman")) +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           panel.border = element_blank(), 
#           axis.line = element_line(colour = "black"),
#           axis.ticks.length=unit(-0.25, "cm"),
#           axis.text.y=element_text(margin=margin(0,15,0,0)), 
#           axis.text.x=element_text(margin=margin(15,0,0,0))) +
#     annotate("text", x = Inf, y = Inf, label = g3label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
#     ylab(expression('Number')) + 
#     xlab(expression('B - R')) + 
#     scale_y_continuous(expand= c(0,0)) 
#     #annotate("text", x = 1.5, y = 0.2, label = galname,size = 6,family="CM Roman", hjust = 0, vjust = 0)

# #plot(mirror.ticks(g))
# gplot <- arrangeGrob(g0, g1, g3, nrow=3, ncol=1)
# ggsave(outfile, gplot, width=8, height=8)
# #ggsave(outfile, g, width=7, height=7)
# #embed_fonts('tmp.pdf', outfile=outfile)

# # circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
# #     r = diameter / 2
# #     tt <- seq(0,2*pi,length.out = npoints)
# #     xx <- center[1] + r * cos(tt)
# #     yy <- center[2] + r * sin(tt)
# #     return(data.frame(x = xx, y = yy))
# # }

# # dat <- circleFun(c(xcen,ycen),extent,npoints = 100)

# # bgccs <- gccdata[gccdata$bmr < m2join$crossing,]
# # rgccs <- gccdata[gccdata$bmr >= m2join$crossing,]

# # s1 <- ggplot(data=bgccs) +
# #       geom_point(aes(x=bx, y=by), color="blue") +
# #       geom_path(data=dat, aes(x,y), linetype=3) +
# #       theme_bw() + 
# #       theme(legend.position="none") + 
# #       theme(text=element_text(family="CM Roman")) +
# #       theme(panel.grid.major = element_blank(), 
# #             panel.grid.minor = element_blank(), 
# #             panel.border = element_blank(), 
# #             axis.line = element_line(colour = "black"),
# #             axis.ticks.length=unit(-0.25, "cm"),
# #             axis.text.y=element_text(margin=margin(0,15,0,0)), 
# #             axis.text.x=element_text(margin=margin(15,0,0,0))) +
# #       # annotate("text", x = Inf, y = Inf, label = g4label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
# #       ylab(expression('Y position (pixels)')) + 
# #       xlab(expression('X position (pixels)'))

# # s2 <- ggplot(data=rgccs) +
# #       geom_point(aes(x=bx, y=by), color="red") +
# #       geom_path(data=dat, aes(x,y), linetype=3) +
# #       theme_bw() + 
# #       theme(legend.position="none") + 
# #       theme(text=element_text(family="CM Roman")) +
# #       theme(panel.grid.major = element_blank(), 
# #             panel.grid.minor = element_blank(), 
# #             panel.border = element_blank(), 
# #             axis.line = element_line(colour = "black"),
# #             axis.ticks.length=unit(-0.25, "cm"),
# #             axis.text.y=element_text(margin=margin(0,15,0,0)), 
# #             axis.text.x=element_text(margin=margin(15,0,0,0))) +
# #       # annotate("text", x = Inf, y = Inf, label = g4label, size = 2, family="CM Roman",hjust = 1, vjust = 1) +
# #       ylab(expression('Y position (pixels)')) + 
# #       xlab(expression('X position (pixels)'))

# # splot <- arrangeGrob(s1, s2, s2, s1, nrow=2, ncol=2)
# # ggsave(outfile2, splot, width=8, height=8)
# # x <- bmr

# # N <- length(x)
# # ks.boot <- rep(0,1000)
# # for (i in 1:1000) {
# #   z <- rbinom(N, 1, foo$lambda[1])
# #   x.b <- z*rnorm(N, foo$mu[1], foo$sigma[1]) + (1-z)*rnorm(N, foo$mu[2], foo$sigma[2])
# #   foo.b <- normalmixEM(x.b, maxit=10000, lambda=foo$lambda, mu=foo$mu, sigma=foo$sigma)
# #   ks.boot[i] <- ks.test(x.b, pmnorm, mu=foo.b$mu, sigma=foo.b$sigma, pmix=foo.b$lambda)$statistic
# # }

# # mean(test$statistic <= ks.boot)