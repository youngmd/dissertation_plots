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
library(ADGofTest)

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

ADiters <- 1000

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
  ad.boot
  #mean(test$statistic <= ad.boot)
}
mwgc <- gccdata
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc <- subset(mwgc, V.R != 0)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc$BR <- mwgc$B.V + mwgc$V.R - (mwgc$A_B - mwgc$A_R)
# mwgc <- subset(mwgc, BR < 1.6)

m <- Mclust(mwgc$BR)

df <- data.frame('BIC'=m$BIC[,])
df$x <- row.names(df)

df3 <- gather(df, Variance, BIC, BIC.E, BIC.V)

g1 <- ggplot(data=df3, aes(x=x,y=BIC)) + geom_point(aes(shape=Variance,color=Variance)) + geom_line(aes(group=Variance,color=Variance)) + 
  theme_bw() + 
  theme(legend.position=c(0.3,0.3)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  xlab('Model Components')

be <- mclustBootstrapLRT(mwgc$BR, modelName='E')
#bv <- mclustBootstrapLRT(mwgc$BR, modelName='V')

bestind <- length(be$G) - 1
secind <- bestind - 1
bestboot <- be$boot[,bestind]
secboot <- be$boot[,secind]
binwa <- 8 * IQR(bestboot) / length(bestboot)^(1/3)
binwb <- 8 * IQR(secboot) / length(secboot)^(1/3)
s <- data.frame('BLRT'=bestboot)
s2 <- data.frame('BLRT'=secboot)

g2a <- ggplot(data=s2, aes(x=BLRT)) + 
  geom_histogram(binwidth=binwb, fill='grey80') + geom_vline(xintercept=be$obs[secind], color='red4', linetype=2) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("Likelihood Ratio (2 vs 3)")

g2b <- ggplot(data=s, aes(x=BLRT)) + 
  geom_histogram(binwidth=binwa, fill='grey80') + geom_vline(xintercept=be$obs[bestind], color='red4', linetype=2) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("Likelihood Ratio (3 vs 4)")

g2 <- arrangeGrob(g2a, g2b, nrow=2, ncol=1)
x <- mwgc$BR
binw <- 1.2 * IQR(x) / length(x)^(1/3)
# binw <- 0.05
m <- Mclust(mwgc$BR, G=3, modelName='E')
dist <- data.frame( BmR = seq(from = min(x)*0.95, to = max(x)*1.05, length.out = 500))
dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sqrt(m$parameters$variance$sigmasq)) * (m$parameters$pro[1] * length(x) * binw)
dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sqrt(m$parameters$variance$sigmasq)) * (m$parameters$pro[2] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sqrt(m$parameters$variance$sigmasq)) * (m$parameters$pro[3] * length(x) * binw)

g3 <- ggplot() + 
  geom_histogram(data=mwgc, aes(x=BR), fill="grey", binwidth=binw) +
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
  annotate("text", x = 1.5, y = 0.9*max(dist$B,dist$R,dist$V), label=galname, size =4, hjust=0, vjust=0.5) +
  scale_y_continuous(expand= c(0,0)) 

m3 <- normalmixEM(mwgc$BR, mu=m$parameters$mean, maxit=10000, lambda=m$parameters$pro, sigma=sqrt(m$parameters$variance$sigmasq), k=3, mean.constr=c(NA, NA, NA), arbmean=TRUE, arbvar=FALSE)
ad.score <- ADGofTest::ad.test(mwgc$BR, pnormmix, m3)
adboot <- bootstrapAD3(mwgc$BR, m3, arbvar=FALSE)

adp <- mean(ad.score$statistic <= adboot)
binw <- 2 * IQR(adboot) / length(adboot)^(1/3)
adb <- data.frame('boot'=adboot)
g4 <- ggplot(data=adb, aes(x=boot)) + 
  geom_histogram(binwidth=binw, fill="#00BFC4") + geom_vline(xintercept=ad.score$statistic, color='#F8766D', linetype=2) +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("Anderson-Darling Score") +
  ylab("Frequency") +
  annotate("text", x = 0.4, y = 0.07*nrow(adb), label=paste0("Bootstrap A-D Test\n\n3 gaussian mixture (N=",nrow(gccdata),")\nEqual variance\np-value=",adp), size = 2, hjust=0, vjust=0.5)



# df <- data.frame(BmR=m$data[,],uncertainty=m$uncertainty, class=m$classification)
# g4 <- ggplot(df, aes(x=BmR, y=uncertainty)) +
#     geom_segment(aes(xend=x, yend=0, color=class)) +
#     theme_bw() + 
#     theme(legend.position="none") + 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           axis.line = element_line(colour = "black"),
#           axis.ticks.length=unit(-0.25, "cm"),
#           axis.text.y=element_text(margin=margin(0,15,0,0)), 
#           axis.text.x=element_text(margin=margin(15,0,0,0))) +
#     scale_y_continuous(limits=c(NA,1), expand=c(0,0))

gplot <- arrangeGrob(g1, g2, g3, g4, nrow=2, ncol=2)
ggsave(datapath, gplot, width=8, height=8)
