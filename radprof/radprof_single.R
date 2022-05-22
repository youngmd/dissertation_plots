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
# Read in the data

cnames <- c("r","N","sd","sderr","lsd","lsderrhi","lsderrlo","good_frac","r_d","NGC","Ncontam","contam_frac")

rp = read.table(paste(datapath,'/corrected_profile.all.out',sep=''), col.names=cnames)
rp <- subset(rp, r < max(outer, 1.5*extent))
rp$Type <- label1

rp
rp$dr <- rp$r^0.25
rp$pr <- log10(rp$r)
rp$lsd <- log10(rp$sd)
rp$lsderrhi <- log10(rp$sd + rp$sderr) - rp$lsd
rp$lsderrlo <- rp$lsd - log10(rp$sd - rp$sderr)
rp_clean = subset(rp, !is.na(lsd))
rp_clean = subset(rp_clean, r < 1.1*extent)

g <- list()

fitlines <- data.frame( r = seq(from = min(rp$r), 
            to = max(rp$r), 
            length.out = 500))

distance = 10^((distmod + 5.) / 5.)
distance_err = 10^((distmod+distmod_err + 5.) / 5.) - distance
extent_kpc = (distance * extent / 3438) / 1000
extent_kpc_err = ((distance+distance_err) * ((extent+extent_err) / 3438) / 1000) - extent_kpc
extent_str <- sprintf("%0.1f ± %0.1f kpc", extent_kpc, extent_kpc_err)

pf <- summary(lm(lsd ~ pr, data=rp_clean, weight=1.0 / lsderrhi^2))
pf$a0 <- pf$coefficients[1,1]
pf$a1 <- pf$coefficients[2,1]

eq <- substitute(italic(slope) == b~"log10(r)+"~a, 
         list(a = format(pf$a0, digits = 2), 
              b = format(pf$a1, digits = 2)))
pow_eq <- as.character(as.expression(eq));    

df <- summary(lm(lsd ~ dr, data=rp_clean, weight=1.0 / lsderrhi^2))
df$a0 <- df$coefficients[1,1]
df$a1 <- df$coefficients[2,1]

dev_eq <- data.frame(V1 = list('tmp'))
eq <- substitute(italic(log10(sigma)) == b~italic(r^{1/4}~"+"~a), 
         list(a = format(df$a0, digits = 3), 
              b = format(df$a1, digits = 3)))
dev_eq$V1 <- as.character(as.expression(eq)); 

dev_label <- as.character(as.expression('R^{1/4}~Fit'))

fitlines$pow <- pf$a0 + pf$a1 * log10(fitlines$r)
fitlines$dev <- df$a0 + df$a1 * fitlines$r^(0.25)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

# mb <- as.data.frame(cbind(x=log10(radprof$r), predict(b, interval='confidence')))

upper = max(rp$sd+rp$sderr, 10^max(fitlines$pow))
outer = max(rp$r)

x = c(0.725*outer, 0.725*outer)
y = c(0.64*upper, 0.56*upper)
type = c(label1, label2)

Types <- data.frame(x, y, type)

rp$Type <- factor(rp$Type, levels = c(label1, label2))
g1 <- ggplot(data=rp) + 
    geom_point(aes(x=r, y=sd, shape=Type, color=Type)) +  
    geom_errorbar(aes(x=r, ymin=sd-sderr,ymax=sd+sderr), width=0.01*outer, size=0.2) +
    geom_vline(aes(xintercept=extent), linetype=4, size=0.2, color='red4') +
    geom_rect(xmin=0.65*outer, xmax=0.98*outer, ymin=0.68*upper, ymax=0.98*upper, fill='white', color='black', alpha=0.95) +
    # geom_point(data=Types, aes(x=x, y=y, shape=type, color=type)) +
    geom_path(data = fitlines, aes( r, 10^pow), linetype=2, size=0.4, color='blue4') +
    geom_path(data = fitlines, aes( r, 10^dev), linetype=3, size=0.4, color='green4') +
    geom_segment(aes(x = 0.7*outer, xend= 0.75*outer, y=0.8*upper, yend=0.8*upper), linetype=2, size=0.4,  color='blue4') +
    geom_segment(aes(x = 0.7*outer, xend= 0.75*outer, y=0.72*upper, yend=0.72*upper), linetype=3, size=0.4, color='green4') +
    geom_blank(aes(x=0-0.05*(extent),y=upper), stat="identity") + 
    geom_blank(aes(x=0-0.05*(extent),y=-0.05*(upper)), stat="identity") + 
    geom_blank(aes(r=1.1*r), stat="identity") +
    geom_hline(aes(yintercept=0), linetype=5, size=0.2) +
    # geom_text(aes(x=extent, label=extent_kpc, family="CM Roman", size=2, y=0.4*max(sd+sderr), angle=90)) +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab(expression('Number arcmin'^-2)) + 
    xlab(expression('Projected radial distance (arcmin)')) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    annotate("text", x = 0.98*extent, y = 0.1*upper, label = extent_str, size = 3, hjust=0, vjust=0.5, family="CM Roman", angle=90) +
    annotate("text", x = 0.95*outer, y = 0.88*upper, label = galname, size = 5, family="CM Roman", hjust = 1, vjust = 0) +
    annotate("text", x = 0.77*outer, y = 0.8*upper, label = 'Power Law Fit', size = 3, family="CM Roman", hjust = 0, vjust = 0.3) +
    annotate("text", x = 0.77*outer, y = 0.72*upper, label = dev_label, parse=TRUE, size = 3, family="CM Roman", hjust = 0, vjust = 0.3) 

#plot(mirror.ticks(g))

g2 <- ggplot(data=rp) + 
    geom_point(aes(x=r, y=sd, shape=Type, color=Type)) + 
    geom_errorbar(aes(x=r, ymin=sd-sderr,ymax=sd+sderr), width=0.02*log10(outer), size=0.2) +
    geom_path(data = fitlines, aes( r, 10^pow), linetype=2, size=0.4, color='blue4') +
    geom_path(data = fitlines, aes( r, 10^dev), linetype=3, size=0.4, color='green4') +
    # geom_text(data=dev_eq,aes(x = 1, y = 1.6,label=V1), parse = TRUE, inherit.aes=FALSE, hjust=0) +
    geom_vline(aes(xintercept=extent), linetype=4, size=0.2, color='red4') +
    geom_hline(aes(yintercept=0), linetype=5, size=0.2) +
    geom_blank(aes(y=1.1*(sd+sderr)), stat="identity") + 
    geom_blank(aes(r=1.1*r), stat="identity") + 
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab(expression(paste('Number arcmin'^-2,''))) + 
    xlab(expression('Projected radial distance (arcmin)')) + 
    scale_y_continuous(trans='asinh',breaks=c(0,1,10,100),labels=c("0",expression(10^0),expression(10^1),expression(10^2)))+
    scale_x_log10(breaks=c(1,2,5,10,20,50,100))

gplot <- arrangeGrob(g1, g2, nrow=2, ncol=1)
# CairoPNG(outfile,height=1024,width=1024)
ggsave(outfile, gplot, width=7, height=6)

pf

df

outer <- 0
# dev.off()
#ggsave(outfile, g, width=7, height=3)
#embed_fonts('tmp.pdf', outfile=outfile)
