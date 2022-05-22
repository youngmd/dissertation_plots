z1 <- read.table('z1sol.dat', header=TRUE)
z2 <- read.table('z05sol.dat', header=TRUE)
z3 <- read.table('z025sol.dat', header=TRUE)
z4 <- read.table('z01sol.dat', header=TRUE)
#z5 <- read.table('z03sol.dat', header=TRUE)
z6 <- read.table('z003sol.dat', header=TRUE)

z1$Z <- "1Z"
z2$Z <- "0.5Z"
z3$Z <- "0.25Z"
z4$Z <- "0.1Z"
#z5$Z <- "0.3Z"
z6$Z <- "0.03Z"

z <- rbind(z1, z2, z3, z4, z6)
z$BR <- z$B - z$R

z$BR <- z$B - z$R
z$BV <- z$B - z$V
library(ggplot2)
library(Cairo)

legend_title <- "SSP"
# 
plotlabels = expression('0.03Z'['\u0298'],
                        '0.1Z'['\u0298'],
                        '0.25Z'['\u0298'],
                        #'0.3Z'['\u0298'],
                        '0.5Z'['\u0298'],
                        Z['\u0298'])

g <- ggplot(z) + 
  geom_path(aes(x=BR, y=Age, color=Z, linetype=Z), size=0.4) + 
  geom_vline(aes(xintercept=1.27), size=0.2, linetype=1, color='#00BA38') + 
  scale_linetype_manual(legend_title, values=c(1,2,5,6,4),
                        labels=plotlabels) +
  scale_color_manual(legend_title, values=c('#00B0F6','#a6611a','#dfc27d','#80cdc1','#018571'), 
                     labels=plotlabels) +
  theme_bw() +
  theme(legend.position=c(0.2,0.8)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=12),
        legend.key = element_blank(),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=10), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=10)) +
  ylab(expression('Age (yrs)')) + 
  xlab(expression('B - R'))  +
  scale_y_log10() +
  annotation_logticks(sides = "lr")

ggsave('isometal.pdf', g, width=6, height=6, device=cairo_pdf)