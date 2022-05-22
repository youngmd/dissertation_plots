library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(extrafont)
library(gtable)
library(gridExtra)
library(igraph)
require(plyr)

insert_minor <- function(major_labs, n_minor) {labs <- 
                              c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
                              labs[1:(length(labs)-n_minor)]}

GCS = read.csv(file="../mastertable/mastertable_gcs.csv",header=TRUE,dec=".",sep=",")


#GCS$Source = factor(GCS$oursurvey, levels=c('This Study','Our survey','Other studies','Galaxy','M31'))

g1 <- ggplot(data=GCS, aes(x=M_V, y=S_N, color=Type, shape=Type)) +
    geom_errorbar(aes(x=M_V, ymin=S_N-S_N_err,ymax=S_N+S_N_err), width=0.05, size=0.2) +
    geom_errorbarh(aes(xmin=M_V-M_V_err,xmax=M_V+M_V_err), height=0.2, size=0.2) +
    geom_point(size=1.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position=c(0.2,0.8),
          legend.key.size=unit(1,"cm"),
          legend.key.height=unit(1,"line"),
          legend.key = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    scale_x_reverse() +
    ylab(expression(paste('Specific Frequency',(S[N])))) + 
    xlab(expression(M[V]^T))

g2 <- ggplot(data=GCS, aes(x=mass, y=T, color=Type, shape=Type)) +
    geom_errorbar(aes(ymin=T-T_err,ymax=T+T_err), width=0.03, size=0.2) +
    geom_errorbarh(aes(xmin=mass-mass_err,xmax=mass+mass_err), height=0.2, size=0.2) +
    geom_point(size=1.5) +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab('Specific Frequency (T)') + 
    xlab(expression(paste('log(Mass/',M[sun],') of Host Galaxy')))

gplot <- arrangeGrob(g1, g2, nrow=2, ncol=1)
ggsave('specfreq.pdf', gplot, width=7, height=7)


GCS <- subset(GCS, !is.na(Extent))
GCS <- subset(GCS, !is.na(mass))
g3 <- ggplot(data=GCS, aes(x=mass, y=Extent)) +
    geom_errorbar(aes(color=Type, ymin=Extent-Extent_err,ymax=Extent+Extent_err), width=0.02, size=0.2) +
    geom_errorbarh(aes(color=Type, xmin=mass-mass_err,xmax=mass+mass_err), height=2, size=0.2) +
    geom_point(aes(color=Type, shape=Type), size=1.5, fill=NA) +
    geom_smooth(method='lm', fullrange=TRUE, formula=y ~ poly(x, 2, raw=TRUE), level=0.99, linetype=1, size=0.5) +
    theme_bw() +
    theme(legend.title=element_blank(),
          legend.position=c(0.2,0.8),
          legend.key.size=unit(1,"cm"),
          legend.key.height=unit(1,"line"),
          legend.key = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab('Spatial Extent of the GC System (kpc)') + 
    xlab(expression(paste('log(Mass/',M[sun],') of Host Galaxy'))) +
    scale_x_continuous(expand=c(0,0))
    

ggsave('extent.pdf', g3, width=7, height=7)
df <- lm(data=GCS, formula=Extent ~ poly(mass, 2, raw=TRUE))
summary(df)
ggsave('extent.pdf', g3, width=7, height=7)