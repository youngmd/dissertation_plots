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

# mwgc <-read.table('harris_clean.txt',header=TRUE)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc <- subset(mwgc, V.R != 0)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc$BR <- mwgc$B.V + mwgc$V.R - (mwgc$A_B - mwgc$A_R)
# mwgc <- subset(mwgc, BR < 1.6)

m <- Mclust(d$bmr)

df <- data.frame('BIC'=m$BIC[,])
df$x <- row.names(df)

df3 <- gather(df, Variance, BIC, BIC.E, BIC.V)

g1 <- ggplot(data=df3, aes(x=x,y=BIC)) + geom_point(aes(shape=Variance,color=Variance)) + geom_line(aes(group=Variance,color=Variance)) + 
    theme_bw() + 
    theme(legend.position=c(0.2,0.3)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    xlab('Model Components')

if(modelName == 'E'){
  be <- mclustBootstrapLRT(bmr, modelName='E')
  blrtlab <- "Equal variance"
} else {
  be <- mclustBootstrapLRT(bmr, modelName='V')
  blrtlab <- "Unequal variance"
}

bestind <- length(be$G)
bestboot <- be$boot[,bestind]
binw <- 8 * IQR(bestboot) / length(bestboot)^(1/3)
s <- data.frame('BLRT'=bestboot)
edge <- max(2*be$obs[bestind],2*sd(s$BLRT))
g2 <- ggplot(data=s, aes(x=bestboot)) + 
    geom_histogram(binwidth=binw, fill='grey80') + geom_vline(xintercept=be$obs[bestind], color='red4', linetype=2) +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(-1,edge)) +
    xlab("Likelihood Ratio") +
    annotate("text", x = 0.5*edge-1, y = 400, label=paste0("Bootstrap Likelihood Ratio Test\n\n",max(be$G)," component vs. ",max(be$G)+1," component\n",blrtlab,"\np-value=",be$p.value[bestind]), size = 3, hjust=0, vjust=0.5)


m <- Mclust(d$bmr, G=3:6, modelName=modelName)

x <- d$bmr
binw <- 1.2 * IQR(x) / length(x)^(1/3)
# binw <- 0.05
dist <- data.frame( BmR = seq(from = min(x)*0.95, to = max(x)*1.05, length.out = 500))
if(modelName == 'E'){
    sd1 = sqrt(m$parameters$variance$sigmasq)
    sd2 = sqrt(m$parameters$variance$sigmasq)
    sd3 = sqrt(m$parameters$variance$sigmasq)
} else {
    sd1 = sqrt(m$parameters$variance$sigmasq[1])
    sd2 = sqrt(m$parameters$variance$sigmasq[2])
    sd3 = sqrt(m$parameters$variance$sigmasq[3])
}

dist$B <- dnorm(dist$BmR,mean=m$parameters$mean[1], sd=sd1) * (m$parameters$pro[1] * length(x) * binw)
dist$V <- dnorm(dist$BmR,mean=m$parameters$mean[2], sd=sd2) * (m$parameters$pro[2] * length(x) * binw)
dist$R <- dnorm(dist$BmR,mean=m$parameters$mean[3], sd=sd3) * (m$parameters$pro[3] * length(x) * binw)


g3 <- ggplot() + 
    geom_histogram(data=d, aes(x=bmr), fill="grey", binwidth=binw) +
    geom_path(data = dist, aes( BmR, B), linetype=2, color='blue4') +
    geom_path(data = dist, aes( BmR, V), linetype=4, color='green4') +
    geom_path(data = dist, aes( BmR, R), linetype=5, color='red4') +
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

df <- data.frame(BmR=m$data[,],uncertainty=m$uncertainty, class=m$classification)
g4 <- ggplot(df, aes(x=BmR, y=uncertainty)) +
    geom_segment(aes(xend=x, yend=0, color=class)) +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    scale_y_continuous(limits=c(NA,1), expand=c(0,0)) +
    xlab(expression('B - R'))

gplot <- arrangeGrob(g1, g2, g3, g4, nrow=2, ncol=2)
ggsave('mm3.pdf', gplot, width=8, height=8)
