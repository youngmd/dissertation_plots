library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
require(plyr)
library(reshape)
library(fontcm)

completeB <- read.table('~/gcc/n5846/gcfinder2/b.complete')
completeV <- read.table('~/gcc/n5846/gcfinder2/v.complete')
completeR <- read.table('~/gcc/n5846/gcfinder2/r.complete')
source('~/gcc/n5846/gcfinder2/photsoln.kvp')

apb = -0.181
apv = -0.215
apr = -0.188
airmass_b <- 1.399
airmass_v <- 1.202
airmass_r <- 1.293
kb = 0.25
kv = 0.15
kr = 0.10

completeB$filt <- 'B'
completeV$filt <- 'V'
completeR$filt <- 'R' 

completeB$mag_o <- completeB$V1 + apb - airmass_b * kb
completeV$mag_o <- completeV$V1 + apv - airmass_v * kv
completeR$mag_o <- completeR$V1 + apr - airmass_r * kr
completeB$mag <- completeB$mag_o + zp_bv + zp_v
completeV$mag <- completeV$mag_o + zp_v
completeR$mag <- completeR$mag_o - zp_vr + zp_v

complete <- rbind(completeB, completeV)
complete <- rbind(complete, completeR)
complete$Filter_s = factor(complete$filt, levels=c('B','V','R'))


g <- ggplot(data=complete, aes(mag, V3, color=Filter_s, shape=Filter_s)) + 
    geom_path(aes(linetype=Filter_s)) +
    #geom_smooth(se=FALSE, aes(linetype = Filter_s), size=0.4) +
    #geom_point(data=completeV, aes(x=mag, y=V3), size=2, shape=21, alpha=0.7, color="green") +
    #geom_smooth(data=completeV, aes(x=V1, y=V2), method=lm, formula=y ~ poly(x,8), linetype=2, color="green") +
    #geom_point(data=completeR, aes(x=mag, y=V3), size=2, shape=20, alpha=0.7, color="red") +
    #geom_smooth(data=completeR, aes(x=V1, y=V2), method=lm, formula=y ~ poly(x,8), linetype=3, color="red") +
    theme_bw() + 
    theme(text=element_text(family="CM Roman")) +
    theme(legend.title=element_blank(),
          legend.position=c(0.2,0.15),
          legend.key.size=unit(1,"cm"),
          legend.key = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    # annotate("text", x = Inf, y = Inf, label = galname, size = 5, family="CM Roman",hjust = 1, vjust = 1) +
    scale_color_manual(values=c("blue3", "green3",
                              "red3")) +
    ylab("Fraction Complete") + 
    xlab("Magnitude") +
    scale_y_continuous(limits=c(0,1.0)) +
    scale_x_continuous(expand = c(0,0))


ggsave("completeness_n5846.pdf", g, width=8, height=8)