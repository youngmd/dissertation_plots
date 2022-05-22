library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)


plotGCLF <- function(galname, datapath, outfile, mean, cutoff){
    source('gclf.R',print.eval=TRUE)
    return(g)
}

cutoff <- 0.7
galname <- "NGC 5846"
datapath <- "~/gcc/n5846/gcfinder3/"
outfile <- "n5846.gclf.pdf"
mean <- 24.58
g1 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

galname <- "NGC 4649"
datapath <- "~/gcc/n4649/gcfinder3.n4649/"
outfile <- "n4649.gclf.pdf"
mean <- 23.73
bw <- 0.4
g2 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

galname <- "NGC 4621"
datapath <- "~/gcc/n4649/gcfinder3.n4621/"
outfile <- "n4621.gclf.pdf"
mean <- 23.91
bw <- 0.5
g3 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

galname <- "NGC 4382"
datapath <- "~/gcc/n4382/gcfinder3/"
outfile <- "n4382.gclf.pdf"
mean <- 23.93
bw <- 0.4
g4 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

galname <- "NGC 1023"
datapath <- "~/gcc/n1023/"
outfile <- "n1023.gclf.pdf"
mean <- 22.96
g5 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

galname <- "NGC 4013"
datapath <- "~/gcc/n4013/"
mean <- 23.9
outfile <- "n4013.gclf.pdf"
g6 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

cutoff <- 0.5
galname <- "NGC 7332"
datapath <- "~/gcc/n7332/"
outfile <- "n7332.gclf.pdf"
mean <- 24.28
g7 <- plotGCLF(galname, datapath, outfile, mean, cutoff)

galname <- "NGC 7339"
datapath <- "~/gcc/n7339/"
outfile <- "n7339.gclf.pdf"
mean <- 24.42
g8 <- plotGCLF(galname, datapath, outfile, mean, cutoff)


df <- data.frame(x = rnorm(100), y = rnorm(100))
df2 <- data.frame(x = rnorm(1), y = rnorm(1))
g9 <- ggplot(df) +
    geom_segment(aes(x = -1, xend= 0, y=0.6, yend=0.6), linetype=2, size=0.3)+
    geom_segment(aes(x = -1, xend= 0, y=0.3, yend=0.3), linetype=3, size=0.3)+
    geom_segment(aes(x = -1, xend= 0, y=0, yend=0), linetype=4, size=0.3)+
    geom_rect(data=df2, xmin=-1, xmax=0, ymin=-0.4, ymax=-0.2)+
    geom_rect(data=df2, xmin=-1, xmax=0, ymin=-0.7, ymax=-0.5, alpha=0.5)+
    geom_blank(aes(x, y)) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    annotate("text", x = 0.1, y = 0.6, label = "1.2", size = 3, hjust=0, vjust=0.5, family="CM Roman") +
    annotate("text", x = 0.1, y = 0.3, label = "1.3", size = 3, hjust=0, vjust=0.5, family="CM Roman") +
    annotate("text", x = 0.1, y = 0.0, label = "1.4", size = 3, hjust=0, vjust=0.5, family="CM Roman") +
    annotate("text", x = 0.1, y = -0.3, label = "Observed", size = 3, hjust=0, vjust=0.5, family="CM Roman") +
    annotate("text", x = 0.1, y = -0.6, label = "Corrected", size = 3, hjust=0, vjust=0.5, family="CM Roman")


gplot <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, nrow=3, ncol=3)
ggsave('gclf.pdf', gplot, width=8, height=8)

# galname <- "NGC 1172"
# datapath <- "~/gcc/Files_for_GCS_Database/n1172/color_sample90.inner_3.45"
# outfile <- "n1172.gclf.pdf"
# source('gclf.R',print.eval=TRUE)

# galname <- "NGC 4406"
# datapath <- "~/gcc/Files_for_GCS_Database/n4406/gc_cand_calib.after_cut.90_percent_sample"
# outfile <- "n4406.mm.pdf"
# source('color_dist_mm.R',print.eval=TRUE)

# galname <- "NGC 1023"
# datapath <- "~/N1023/"
# outfile <- "n1023.cc.pdf"
# source('colorcolor.R',print.eval=TRUE)

# galname <- "NGC 7332"
# datapath <- "~/gcc/n7332/"
# outfile <- "n7332.cc.pdf"
# source('colorcolor.R',print.eval=TRUE)

# galname <- "NGC 4013"
# datapath <- "~/gcc/n4013/"
# outfile <- "n4013.cc.pdf"
# source('colorcolor.R',print.eval=TRUE)