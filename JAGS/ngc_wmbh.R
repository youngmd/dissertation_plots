library(rjags)
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
require(runjags)

source('plotfuncs.R', print.eval=TRUE)

GCS = read.csv(file="../gcs2.csv",header=TRUE,dec=".",sep=",")
GCS = subset(GCS, !is.na(lg_M_B)) 
GCS = subset(GCS, !is.na(lg_M_B_err_high)) 
GCS$lg_M_B_err <- GCS$lg_M_B_err_high
GCS = subset(GCS, !is.na(Mdyn)) 
GCS = subset(GCS, !is.na(N_GC)) 
GCS = subset(GCS, !is.na(Mdyn_err))
GCS = subset(GCS, !is.na(N_GC_err)) 
GCS = subset(GCS, !is.na(M_V)) 
GCS = subset(GCS, !is.na(M_V_err))
GCS = subset(GCS, !is.na(M_K)) 
GCS = subset(GCS, !is.na(sig_e))
GCS = subset(GCS, !is.na(R_e))

burnin <- 20000
sample <- 30000
thin <- 3
# Allocate 3 cores in the computer
#
##
datapath = "FiguresMBH/"

calcD <- function(Dstats, pred, D1, D2){
  DIC <- sum(D2$deviance) + sum(D2$penalty)
  Dstats2 <- data.frame('pred'=pred,'dispersion'=D1, 'DIC'=DIC)
  Dstats <- rbind(Dstats, Dstats2)
  return(Dstats)
}

Dstats <- data.frame('pred'=character(),'dispersion'=double(), 'DIC'=double())

ylabel <- expression(N[GC])
N_low<-GCS$N_GC-GCS$N_GC_err
N_low[N_low<0]<-0
GCS$Type <- factor(GCS$group_type, levels = c("E", "Spiral", "S0"))

GCS$L_V <- (4.83 - GCS$M_V)/2.5
GCS$L_V_err <- (4.83 - GCS$M_V + GCS$M_V_err)/2.5 - GCS$L_V

GCS$L_K <- (3.27 - GCS$M_K)/2.5
GCS$L_K_err <- (3.27 - GCS$M_K + GCS$M_K_err)/2.5 - GCS$L_K


df <- data.frame(y = GCS$N_GC, y_err = GCS$N_GC_err)

summaries <- list()
outputs <- list()

preds <- c("lg_M_B","Mdyn","mass","R_e","L_K","L_V")

for (p in preds) {
  pind <- which(preds==p)
  gc()
  df$x <- GCS[[p]]
  p_err <- paste(p, "_err", sep="")
  df$up_x <- GCS[[p_err]]
  print(paste("Working on",p))
  cl <- makeCluster(3)
  source('nb.R',print.eval=TRUE)
  model = p
  summaries[[model]] = summary
  outputs[[model]] = pred.NB2errx
  stopCluster(cl)
  gc()
  # DIC = summary$dic$dic
  # bod2[counter,]$DIC = DIC
  Dstats <- calcD(Dstats, model, 0, dicsamples.nb)
  
  rm(jags.neg)
  rm(dicsamples.nb)
  rm(jagssamples.nb)
  rm(jags.DIC)
  # # rm(Pres)
  # rm(pred.NBerrx)
  rm(pred.NBerrx2)
}

save(summaries, file='MBHsummaries_autorun.Robj')
save(outputs, file='MBHoutputs_autorun.Robj')
save(Dstats, file='MBHdstats_autorun.Robj')
# source('Sig_e.R', print.eval=TRUE)
# xlabel <- expression(sigma~~(km/s))
# GCS$x <- GCS$sig_e
# GCS$xlow <- GCS$sig_e_err
# GCS$xhigh <- GCS$sig_e_err
# pred.NB2errx$x <- pred.NB2errx$sigex
# g1 <- plotfit(GCS, pred, xlabel, ylabel)
# ggsave(paste(datapath,'Sig.pdf'), g1, width=8, height=8, device=cairo_pdf)
# Dstats <- calcD(Dstats, 'Sig', Dispersion, dicsamples.nb)

# source('MV.R', print.eval=TRUE)
# xlabel <- expression(M[V])
# GCS$x <- GCS$M_V
# GCS$xlow <- GCS$M_V_err
# GCS$xhigh <- GCS$M_V_err
# pred.NB2errx$x <- pred.NB2errx$MV_Tx
# g2 <- plotfit(GCS, pred, xlabel, ylabel)
# g2 <- g2 + scale_x_reverse(expand=c(0,0))
# ggsave(paste(datapath,'MV.pdf'), g2, width=8, height=8, device=cairo_pdf)
# Dstats <- calcD(Dstats, 'MV', Dispersion, dicsamples.nb)

# source('R_e.R', print.eval=TRUE)
# xlabel <- expression(log~R[e]~~(kpc))
# GCS$x <- GCS$R_e
# GCS$xlow <- GCS$R_e_err
# GCS$xhigh <- GCS$R_e_err
# pred.NB2errx$x <- pred.NB2errx$Rex
# g3 <- plotfit(GCS, pred, xlabel, ylabel)
# ggsave(paste(datapath,'Re.pdf'), g3, width=8, height=8, device=cairo_pdf)
# Dstats <- calcD(Dstats, 'R_e', Dispersion, dicsamples.nb)

# source('Mass.R', print.eval=TRUE)
# xlabel <- expression(log~M[stellar]/M['\u0298'])
# GCS$x <- GCS$mass
# GCS$xlow <- GCS$mass_err
# GCS$xhigh <- GCS$mass_err
# pred.NB2errx$x <- pred.NB2errx$massx
# g4 <- plotfit(GCS, pred, xlabel, ylabel)
# ggsave(paste(datapath,'Mass_stellar.pdf'), g4, width=8, height=8, device=cairo_pdf)
# # Dstats <- calcD(Dstats, 'M_stellar', Dispersion, dicsamples.nb)
# df <- data.frame(y = GCS$N_GC, y_err = GCS$N_GC_err)
# # df$y <- GCS$N_GC
# # df$y_err <- GCS$N_GC_err
# df$x <- GCS$Mdyn
# df$downx <- GCS$Mdyn_err
# df$up_x <- GCS$Mdyn_err
# df$z <- GCS$mass
# df$downz <- GCS$mass_err
# df$up_z <- GCS$mass_err
# source('nb2.R', print.eval=TRUE)
# #xlabel <- expression(log~M[dyn]/M['\u0298'])
# GCS$x <- GCS$Mdyn
# GCS$xlow <- GCS$Mdyn_err
# GCS$xhigh <- GCS$Mdyn_err
# #pred.NB2errx$x <- pred.NB2errx$dfx
# #g5 <- plotfit(df, pred, xlabel, ylabel)
# #ggsave(paste(datapath,'Mdyn2.pdf'), g5, width=8, height=8, device=cairo_pdf)
# Dstats <- calcD(Dstats, 'M_dyn', Dispersion, dicsamples.nb)

# source('MBH.R', print.eval=TRUE)
# xlabel <- expression(log~M[SMBH]/M['\u0298'])
# GCS$x <- GCS$lg_M_B
# GCS$xlow <- GCS$lg_M_B_err_low
# GCS$xhigh <- GCS$lg_M_B_err_high
# pred.NB2errx$x <- pred.NB2errx$MBHx
# g6 <- plotfit(GCS, pred, xlabel, ylabel)
# ggsave(paste(datapath,'MBH.pdf'), g6, width=8, height=8, device=cairo_pdf)
# Dstats <- calcD(Dstats, 'MBH', Dispersion, dicsamples.nb)


# E <- subset(GCS, Type == 'E')
# Spiral <- subset(GCS, Type == 'Spiral')
# S0 <- subset(GCS, Type == 'S0')

# N_E <- nrow(E)
# N_Spiral <- nrow(Spiral)
# N_S0 <- nrow(S0)

# E_label = paste("E (N=",N_E,")", sep = "")
# Spiral_label = paste("Spiral (N=",N_Spiral,")", sep = "")
# S0_label = paste("S0 (N=",N_S0,")", sep = "")
# df <- data.frame(x = rnorm(100), y = rnorm(100))
# df2 <- data.frame(x = rnorm(1), y = rnorm(1))
# g7 <- ggplot(df) +
#     geom_rect(data=df2, xmin=-0.6, xmax=1.5, ymin=-1.2, ymax=1.0, color="black", size=0.5, fill=NA) +
#     geom_point(aes(x = -0.3, y=0.8), shape=0, size=2, color='#e41a1c')+
#     geom_point(aes(x = -0.3, y=0.5), shape=1, size=2, color='#377eb8')+
#     geom_point(aes(x = -0.3, y=0.2), shape=5, size=2, color='#4daf4a')+
#     geom_rect(data=df2, xmin=-0.5, xmax=-0.1, ymin=-0.5, ymax=-0.3, fill="grey70", alpha=0.75)+
#     geom_rect(data=df2, xmin=-0.5, xmax=-0.1, ymin=-0.8, ymax=-0.6, fill="grey70", alpha=0.50)+
#     geom_rect(data=df2, xmin=-0.5, xmax=-0.1, ymin=-1.1, ymax=-0.9, fill="grey70", alpha=0.25)+
#     geom_segment(aes(x = -0.5, xend=-0.1, y=-0.1, yend=-0.1), linetype=3, size=0.3)+
#     # geom_rect(data=df2, xmin=-1, xmax=0, ymin=-0.4, ymax=-0.2)+
#     # geom_rect(data=df2, xmin=-1, xmax=0, ymin=-0.7, ymax=-0.5, alpha=0.5)+
#     geom_blank(aes(x, y)) + 
#     scale_y_continuous(expand = c(0,0)) +
#     scale_x_continuous(expand = c(0,0)) +
#     scale_color_brewer(palette='Set1') +
#     theme_bw() + 
#     theme(axis.line=element_blank(),axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none",
#           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank(),plot.background=element_blank()) +
#     annotate("text", x = 0, y = 0.8, label =E_label, size = 3, hjust=0, vjust=0.5) +
#     annotate("text", x = 0, y = 0.5, label =Spiral_label, size = 3, hjust=0, vjust=0.5) +
#     annotate("text", x = 0, y = 0.2, label =S0_label, size = 3, hjust=0, vjust=0.5) +
#     annotate("text", x = 0, y = -0.1, label ="Expected Value", size = 3, hjust=0, vjust=0.5) +
#     annotate("text", x = 0, y = -0.4, label ="50% pred. interval", size = 3, hjust=0, vjust=0.5) +
#     annotate("text", x = 0, y = -0.7, label ="95% pred. interval", size = 3, hjust=0, vjust=0.5) +
#     annotate("text", x = 0, y = -1.0, label ="99% pred. interval", size = 3, hjust=0, vjust=0.5)



# splot <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, nrow=4, ncol=2)
# ggsave(paste(datapath,'Ngc_wmbh_plots.pdf'), splot, width=8, height=10.5, device=cairo_pdf)

# save.image()

# p <- ggplot(Dstats, aes(x=pred, y=DIC, fill=pred)) + 
#     geom_bar(stat="identity") + 
#     geom_text(aes(x=pred,
#                    y=DIC + 0.3 * sign(DIC),
#                    label=format(Dstats$DIC, digits=4)),
#                hjust=0, 
#                size=3,
#                color=rgb(100,100,100, maxColorValue=255)) +
#     labs(x="",
#             y="",
#             title=expression(DIC~of~N[GC]~Predictors)) +
#     theme(axis.text.x =
#          element_text(size  = 10,
#                       angle = 0,
#                       hjust = 1,
#                       vjust = 1),
#          axis.text.y =
#          element_text(hjust = 0.5),
#          panel.background = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.ticks  = element_blank(),
#          axis.line   = element_line(colour=NA),
#          axis.line.x = element_line(colour="grey80")) +
#     theme(legend.position="none") + 
#     coord_flip() +
#     scale_x_discrete(labels=parse(text=c("sigma","M[V]","R[e]","Mass[stellar]","Mass[dyn]","Mass[SMBH]"))) +
#     scale_fill_economist() 

# CairoPDF(paste(datapath,'Ngc_all_DIC.pdf',sep=''),height=5,width=9)
# p
# dev.off()

# p2 <- ggplot(Dstats, aes(x=pred, y=dispersion - 1, fill=pred)) + 
#     geom_bar(stat="identity", position = "identity") + 
#     geom_text(aes(x=pred,
#                    y=(dispersion - 1) + 0.005 * sign(dispersion - 1 ),
#                    label=format(Dstats$dispersion, digits=3)),
#                hjust=ifelse(Dstats$dispersion > 1,0,1),  
#                size=3,
#                color=rgb(100,100,100, maxColorValue=255)) +
#     labs(x="",
#             y="",
#             title=expression(Dispersion~of~N[GC]~Predictors)) +
#     theme(axis.text.x =
#          element_text(size  = 10,
#                       angle = 0,
#                       hjust = 1,
#                       vjust = 1),
#          axis.text.y =
#          element_text(hjust = 0.5),
#          panel.background = element_blank(),
#          panel.grid.minor = element_blank(),
#          axis.ticks  = element_blank(),
#          axis.line   = element_line(colour=NA),
#          axis.line.x = element_line(colour="grey80")) +
#     theme(legend.position="none") + 
#     coord_flip() +
#     scale_x_discrete(labels=parse(text=c("sigma","M[V]","R[e]","Mass[stellar]","Mass[dyn]","Mass[SMBH]"))) +
#     scale_fill_economist() + 
#     scale_y_continuous(breaks=seq(-1,1,0.5), labels=seq(0,2,0.5), limits=c(-0.5,0.5))

# CairoPDF(paste(datapath,'Ngc_all_disp.pdf',sep=''),height=5,width=9)
# p2
# dev.off()

