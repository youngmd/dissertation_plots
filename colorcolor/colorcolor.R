library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(extrafont)
library(gtable)
# Read in the data


#from here http://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
  {
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
      } else {
        x[1:nth != 1]
        }
    } else {
      if(empty) {
        x[1:nth != 1] <- ""
        x
        } else {
          x[1:nth == 1]
        }
    }
}

cnames <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr","galcenrad")

galdata = read.table('new_redshifted_colors.out', header=TRUE)
ptdata = read.table(paste(datapath,'/gc_cand_calib.dat',sep=''))
gccdata = read.table(paste(datapath,'/keep.out',sep=''), col.names=cnames)

box_x <- c(0.5121,0.5879,1.028,0.9521)
box_y <- c(0.4998,0.3044,0.4898,0.6852)
box_data = data.frame(x=box_x,y=box_y)

custom_breaks_x <- seq(0, 2.2, 0.1)
custom_breaks_y <- seq(0.0, 1.5, 0.1)

# ps<-read.table(paste(datapath,'/photo_soln-1.txt',sep=''),sep='#',comment.char="",strip.white=TRUE)
# mu_bv <- ps$V1[which(ps$V2=='mu_bv')]
# mu_bv_err <- ps$V1[which(ps$V2=='mu_bv_err')]
# zp_bv_err <- ps$V1[which(ps$V2=='zp_bv_err')]
# mu_vr <- ps$V1[which(ps$V2=='mu_vr')]
# mu_vr_err <- ps$V1[which(ps$V2=='mu_vr_err')]
# zp_vr_err <- ps$V1[which(ps$V2=='zp_vr_err')]
# kb <- ps$V1[which(ps$V2=='kb')]
# kv <- ps$V1[which(ps$V2=='kv')]
# kr <- ps$V1[which(ps$V2=='kr')]

# gccdata$bo <- gccdata$bi - kb * gccdata$amb
# gccdata$vo <- gccdata$vi - kv * gccdata$amv
# gccdata$ro <- gccdata$ri - kr * gccdata$amr

# gccdata$bmverr <- ((mu_bv_err^2)/(mu_bv^2)) + (gccdata$bi_err^2 + gccdata$vi_err^2)/((gccdata$bo-gccdata$vo)^2)
# gccdata$bmverr <- gccdata$bmverr*(mu_bv*(gccdata$bo-gccdata$vo)^2)
# gccdata$bmverr <- gccdata$bmverr + zp_bv_err^2
# gccdata$bmverr <- sqrt(gccdata$bmverr)
# gccdata$vmrerr <- ((mu_vr_err^2)/(mu_vr^2)) + (gccdata$vi_err^2 + gccdata$ri_err^2)/((gccdata$vo-gccdata$ro)^2)
# gccdata$vmrerr  <- gccdata$vmrerr*(mu_vr*(gccdata$vo-gccdata$ro)^2)
# gccdata$vmrerr  <- gccdata$vmrerr + zp_vr_err^2
# gccdata$vmrerr <- sqrt(gccdata$vmrerr)

gccdata$bmverr = sqrt(gccdata$bi_err^2 + gccdata$vi_err^2)
gccdata$vmrerr = sqrt(gccdata$vi_err^2 + gccdata$ri_err^2)

x_err <- mean(gccdata$bmverr)
y_err <- mean(gccdata$vmrerr)

sample_err = data.frame(x=2.0,y=0.2, ymin=0.2-y_err, ymax=0.2+y_err, xmin=2.0-x_err, xmax=2.0+x_err)

sample_err
# bmverr(i) = ((mu_bv_err**2.)/(mu_bv**2.)) + 
#      &((bi_err(i)**2. + vi_err(i)**2.)/((bo(i)-vo(i))**2.))
#          bmverr(i) = bmverr(i)*((mu_bv*(bo(i)-vo(i)))**2.)
#          bmverr(i) = bmverr(i) + zp_bv_err**2.

# vmrerr(i) = ((mu_vr_err**2.)/(mu_vr**2.)) + 
#      &((vi_err(i)**2. + ri_err(i)**2.)/((vo(i)-ro(i))**2.))
#          vmrerr(i) = vmrerr(i) * ((mu_vr*(vo(i)-ro(i)))**2.)
#          vmrerr(i) = vmrerr(i) + zp_vr_err**2.

# create the plot
g <- ggplot() + 
    geom_point(data = galdata, aes( B.V, V.R, color=factor(galdata$Type))) +
    geom_path(data = galdata, aes( B.V, V.R, color=factor(galdata$Type))) +
    geom_point(data = ptdata, aes(V18, V19), shape=22, size=1, alpha=0.5) + 
    geom_point(data = gccdata, aes(bmv, vmr), size=2) +
    geom_polygon(box_data, mapping=aes(x=x,y=y), fill=NA, color='red',size=0.4) +
    geom_point(data = sample_err, aes(x,y), size=2) +
    geom_errorbar(data = sample_err, aes(x=x,y=y,ymin = ymin,ymax = ymax), width=0.03) +
    geom_errorbarh(data = sample_err, aes(x=x,y=y,xmin = xmin,xmax = xmax), height=0.02) +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    xlab(expression('(B-V)'[0])) + 
    ylab(expression('(V-R)'[0])) + 
    scale_x_continuous(limits = c(0, 2.2), expand=c(0,0), breaks=custom_breaks_x, labels= every_nth(custom_breaks_x, 5, inverse=TRUE)) +
    scale_y_continuous(limits = c(0.0, 1.5), expand=c(0,0), breaks=custom_breaks_y, labels= every_nth(custom_breaks_y, 5, inverse=TRUE)) +
    geom_segment(aes(x =  0.2, y = 1.2, xend = 0.5, yend = 1.38),color='red', arrow = arrow(length = unit(0.25, "cm"))) +
    annotate("text", x = 1.45, y = 1.43, label = "E/S0",size = 4,family="CM Roman") + 
    annotate("text", x = 1.22, y = 1.34, label = "Sab",size = 4,family="CM Roman") + 
    annotate("text", x = 0.9, y = 1.19, label = "Sbc",size = 4,family="CM Roman") + 
    annotate("text", x = 0.7, y = 0.97, label = "Scd",size = 4,family="CM Roman") + 
    annotate("text", x = 0.43, y = 0.73, label = "Irr",size = 4,family="CM Roman") +
    annotate("text", x = 0.1, y = 0.35, label = "z = 0", size = 4, family="CM Roman", hjust = 0, vjust = 0) +
    annotate("text", x = 0.1, y = 0.7, label = "z = 0.7", size = 4, family="CM Roman", hjust = 0, vjust = 0) +
    annotate("text", x = 1.5, y = 0.2, label = galname,size = 6,family="CM Roman", hjust = 0, vjust = 0)

#plot(mirror.ticks(g))

ggsave(outfile, g, width=7, height=7)
#embed_fonts('tmp.pdf', outfile=outfile)
