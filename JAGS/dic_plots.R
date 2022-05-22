library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)
library(plyr)
library(MASS)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(gtable)
library(parallel)
library(viridis)



load('Alldstats_autorun.Robj')
DstatsAll <- Dstats
rm(Dstats)
p1 <- ggplot(DstatsAll, aes(x=pred, y=DIC, fill=pred)) +
  geom_bar(stat="identity") +
  geom_text(aes(x=pred, 
                y=DIC + 0.3 * sign(DIC),
                label=format(DstatsAll$DIC, digits=4)),
                hjust=0,
                size=3,
                color=rgb(100,100,100, maxColorValue=255)) +
  labs( x="", 
        y="", 
        title=expression(DIC~of~N[GC]~Predictors)) +
  theme(axis.text.x =
          element_text(size  = 10,
                       angle = 0,
                       hjust = 1,
                       vjust = 1),
        axis.text.y =
          element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks  = element_blank(),
        axis.line   = element_line(colour=NA),
        axis.line.x = element_line(colour="grey80")) +
  theme(legend.position="none") +
  coord_flip() +
  scale_x_discrete(labels=parse(text=c("log~M[dyn]/M['\u0298']","log~M['*']/M['\u0298']","R[eff]","L[K]","L[V]"))) +
  scale_fill_brewer()

# ggsave('Figures/ngc_all_dic.pdf', p, width=8, height=3, device=cairo_pdf)


load('MBHdstats_autorun.Robj')
DstatsBH <- Dstats
rm(Dstats)

p2 <- ggplot(DstatsBH, aes(x=pred, y=DIC, fill=pred)) +
  geom_bar(stat="identity") +
  geom_text(aes(x=pred, 
                y=DIC + 0.3 * sign(DIC),
                label=format(DstatsBH$DIC, digits=4)),
            hjust=0,
            size=3,
            color=rgb(100,100,100, maxColorValue=255)) +
  labs( x="", y="", title = expression(DIC~of~N[GC]~Predictors~~SMBH~sample)) +
  theme(axis.text.x =
          element_text(size  = 10,
                       angle = 0,
                       hjust = 1,
                       vjust = 1),
        axis.text.y =
          element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks  = element_blank(),
        axis.line   = element_line(colour=NA),
        axis.line.x = element_line(colour="grey80")) +
  theme(legend.position="none") +
  coord_flip() +
  scale_fill_brewer() +
  scale_x_discrete(labels=parse(text=c("log~M[BH]/M['\u0298']","log~M[dyn]/M['\u0298']","log~M['*']/M['\u0298']","R[eff]","L[K]","L[V]" )))

splot <- arrangeGrob(p1, p2, nrow=2, ncol=1)
ggsave('Figures/ngc_dic.pdf', p, width=8, height=8, device=cairo_pdf)


