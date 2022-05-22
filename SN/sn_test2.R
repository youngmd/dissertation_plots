library(ggplot2)

GCS = read.csv(file="../mastertable/mastertable_gcs.csv",header=TRUE,dec=".",sep=",")

GCS <- data.frame(M_V = GCS$M_V, S_N = GCS$S_N, 
                  M_V_err = GCS$M_V_err, S_N_err = GCS$S_N_err,
                  mass = GCS$mass, mass_err = GCS$mass_err,
                  T_N = GCS$T, T_N_err = GCS$T_err, 
                  Extent = GCS$Extent, Extent_err= GCS$Extent_err)

# g5 <- ggplot(data=GCS) +
#   geom_point(aes(x=M_V, y=S_N), size=2) +
#   geom_errorbar(aes(ymin=S_N-S_N_err,ymax=S_N+S_N_err), width=0.1, size=0.2) +
#   geom_errorbarh(aes(xmin=M_V-M_V_err,xmax=M_V+M_V_err), height=0.1, size=0.2) +
#   theme_bw() +
#   theme(legend.position="none") + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.ticks.length=unit(-0.25, "cm"),
#         axis.text.y=element_text(margin=margin(0,15,0,0)), 
#         axis.text.x=element_text(margin=margin(15,0,0,0))) +
#   #scale_x_reverse(expand = c(0,0), breaks=seq(-20,-24,-0.2), labels=insert_minor(seq(-20,-24,-1), 4)) +
#   #scale_y_continuous(expand = c(0,0), breaks=seq(0,10,0.5), labels=insert_minor(seq(0,10,1), 1)) + 
#   ylab(expression(paste('Specific Frequency',(S[N])))) + 
#   xlab(expression(M[V]^T))
# 

g2 <- ggplot() +
  geom_point(data=GCS, aes(x=mass, y=T_N), size=2) +
  geom_errorbar(data=GCS, aes(x=mass, y=T_N, ymin=T_N-T_N_err,ymax=T_N+T_N_err), width=0.1, size=0.2) +
  geom_errorbarh(data=GCS, aes(x=T_N, y=T_N, xmin=mass-mass_err, xmax=mass+mass_err), height=0.1, size=0.2) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)),
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  #scale_x_continuous(expand = c(0,0), breaks=seq(10,13,0.1), labels=insert_minor(seq(10,13,0.5), 4)) +
  #scale_y_continuous(expand = c(0,0), breaks=seq(0,10,0.5), labels=insert_minor(seq(0,10,1), 1)) +
  ylab('Specific Frequency (T)') +
  xlab(expression(paste('log(Mass/',M[sun],') of Host Galaxy')))
