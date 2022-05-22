library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(extrafont)
library(gtable)
# Read in the data

cnames <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr","galcenrad")

gccdata = read.table(paste(datapath,'/gc_cand_calib.90.out',sep=''), col.names=cnames)

dist <- data.frame( BmR = seq(from = 0.75, to = 1.75, length.out = 500))

dist$distB <- dnorm(dist$BmR,mean=1.117222, sd=0.1102) * (257 * 0.05)
dist$distR <- dnorm(dist$BmR,mean=1.433297, sd=0.125) * (84 * 0.05)
# dist$dist1.3 <- dnorm(dist$Vmag,mean=23.93, sd=1.3) * 565.871972
# dist$dist1.4 <- dnorm(dist$Vmag,mean=23.93, sd=1.4) * 541.143759

g <- ggplot() + 
    geom_histogram(data=gccdata, aes(x=(bmv + vmr)), binwidth=0.05) +
    # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=0.4, alpha=0.5) +
    geom_path(data = dist, aes( BmR, distB), linetype=2) +
    geom_path(data = dist, aes( BmR, distR), linetype=3) +
    # geom_path(data = dist, aes( Vmag, dist1.3), linetype=3) +
    # geom_path(data = dist, aes( Vmag, dist1.4), linetype=4) +
    # geom_blank(data=dist, aes(y=1.1*dist1.2), stat="identity") + 
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab(expression('Number')) + 
    xlab(expression('B - R')) + 
    scale_y_continuous(expand = c(0,0)) 
    #annotate("text", x = 1.5, y = 0.2, label = galname,size = 6,family="CM Roman", hjust = 0, vjust = 0)

#plot(mirror.ticks(g))

ggsave(outfile, g, width=7, height=7)
#embed_fonts('tmp.pdf', outfile=outfile)
