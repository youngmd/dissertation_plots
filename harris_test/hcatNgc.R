library(ggplot2)
library(tidyr)
library(Cairo)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(plyr)
library(MASS)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

hcat <- read.csv('harris.csv', header=TRUE)

hcat$Mdyn <- log10(930 * (hcat$R_e*1000) * hcat$sig_e^2)

hcat$L_K <- (3.27 - hcat$M_K)/2.5


# hcat$pNGC <- exp(-26.283 + 1.101*hcat$L_K + 1.842*hcat$Mdyn)
# hcat$pNGClow <- 0
# hcat$pNGChigh <- 0
# load('samp.Robj')
# load('mdyn_lk_pred2All.Robj')
# 
# for(i in 1:nrow(hcat)){
#   print(hcat[i,'Galaxy'][1])
#   lk <- hcat[i,'L_K'][1]
#   md <- hcat[i,'Mdyn'][1]
#   mu <-  exp(samp$beta.0 + samp$beta.1*lk + samp$beta.2*md)
#   p <- samp$size/(samp$size+mu)
#   ngcsamps <- dnbinom(mu, size=samp$size, prob=p)
#   hcat[i,'pNGC'] <- mean(ngcsamps)
#   quants <- quantile(ngcsamps, probs=c(0.01,0.99))
#   hcat[i,'pNGClow'] <- quants[1]
#   hcat[i,'pNGChigh'] <- quants[2]
# }

g <- ggplot(data=pred.NB2errx, aes(x=mean, y=NGC)) + 
  geom_point(aes(color=Type, shape=Type), size=1) + 
  geom_errorbar(aes(color=Type, ymin=NGC-NGC_err,ymax=NGC+NGC_err), size=0.2) + 
  geom_errorbarh(aes(color=Type, xmin=lwr1,xmax=upr1), size=0.2) +
  geom_abline(aes(intercept=0, slope=1), size=0.3, linetype=2) + 
  #geom_point(data=pred.NB2err, aes(x=mean, y=NGC), size=1, color='black') +
  #geom_errorbarh(data=pred.NB2err, aes(x=mean, y=NGC, xmin=lwr3,xmax=upr3), size=0.3, color='black') +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=12),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=8),
        axis.text.x=element_text(margin=margin(5,0,0,0), size=8)) +
  xlab(expression(Predicted~N[GC])) +
  ylab(expression(Observed~N[GC])) +
  scale_y_continuous(trans='asinh',breaks=c(0,1,10,100,1000,10000,100000),labels=c("0",expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)), limits=c(0,NA)) +
  scale_x_continuous(trans='asinh',breaks=c(0,1,10,100,1000,10000,100000),labels=c("0",expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)), limits=c(0,NA))
  
print(nrow(subset(pred.NB2errx, NGC < upr1 & NGC > lwr1 & mean < 100)))
print(nrow(subset(pred.NB2errx, NGC < upr2 & NGC > lwr2 & mean < 100)))
print(nrow(subset(pred.NB2errx, NGC < upr3 & NGC > lwr3 & mean < 100)))
print(nrow(subset(pred.NB2errx, mean < 100)))

# g2 <-  ggplot(data=hcat, aes(x=exp(L_K+Mdyn), y=pNGC)) +
#   geom_point(aes(color=Type, shape=Type), size=0.5) +
#   geom_point(aes(y=N_GC, color=Type, shape=Type), size=0.5) +
#   scale_x_log10() +
#   scale_y_log10() +
#   geom_abline(aes(intercept=0, slope=1)) +
#   geom_errorbar(aes(color=Type, ymin=pNGC-pNGC_err,ymax=pNGC+pNGC_err))

