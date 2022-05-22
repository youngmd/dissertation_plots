z1 <- read.table('z1.br125.dat', header=TRUE)
z2 <- read.table('z05.br125.dat', header=TRUE)
z3 <- read.table('z025.br125.dat', header=TRUE)
z4 <- read.table('z01.br125.dat', header=TRUE)

z1$SSP <- "1"
z1$Age <- "2.3Gyr"
z2$SSP <- "0.5Z"
z2$Age <- "2.7Gyr"
z3$SSP <- "0.25"
z3$Age <- "5.2Gyr"
z4$SSP <- "0.10Z"
z4$Age <- "10Gyr"

z <- rbind(z1, z2, z3, z4)

z$BR <- z$B - z$R
z$BV <- z$B - z$V
library(ggplot2)
library(Cairo)

legend_title <- "SSP (B-R=1.27)"

plotlabels = expression('0.1Z'['\u0298']~'@9.9Gyr',
                        '0.25Z'['\u0298']~'@5.2Gyr',
                        '0.5Z'['\u0298']~'@2.7Gyr',
                        Z['\u0298']~'@2.3Gyr')
g <- ggplot(z) + 
  geom_path(aes(x=BR, y=R, color=SSP, linetype=SSP), size=0.4) + 
  scale_y_reverse(limits=c(15,-3)) + 
  scale_x_continuous(limits=c(0.5,3.5)) + 
  scale_linetype_manual(legend_title, values=c(2,5,6,4),
                        labels=plotlabels) +
  scale_color_manual(legend_title, values=c('#a6611a','#dfc27d','#80cdc1','#018571'), 
                     labels=plotlabels) +
  theme_bw() +
  theme(legend.position=c(0.2,0.2)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.08, "cm"),
        axis.title=element_text(size=12),
        legend.key = element_blank(),
        axis.text.y=element_text(margin=margin(0,5,0,0), size=10), 
        axis.text.x=element_text(margin=margin(5,0,0,0), size=10)) +
  ylab(expression('R')) + 
  xlab(expression('B - R')) 

ggsave('ssp125.pdf', g, width=6, height=6, device=cairo_pdf)