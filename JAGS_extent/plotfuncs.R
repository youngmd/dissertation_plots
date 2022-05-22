
asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

plotfit <- function(GCS, pred, xlabel, ylabel){
  range <- abs(max(GCS$x)) - abs(min(GCS$x))
  whisker <- 0.015 * range
  g <- ggplot(GCS,aes(x=x,y=y))+
  geom_ribbon(data=pred.Norm2errx,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.75, fill="grey70") +
  geom_ribbon(data=pred.Norm2errx,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.50, fill="grey70") +
  geom_ribbon(data=pred.Norm2errx,aes(x=x,y=mean,ymin=lwr3, ymax=upr3), alpha=0.25, fill="grey70") +
  geom_point(aes(x=x, y=y, colour=Type,shape=Type),size=2,alpha=0.8)+
  geom_errorbar(aes(colour=Type,ymin=y-ylow,ymax=y+yhigh),alpha=0.6,width=whisker)+
  geom_errorbarh(aes(colour=Type,xmin=x-xlow, xmax=x+xhigh),alpha=0.6,height=0.15)+
  geom_line(data=pred.Norm2errx,aes(x=x,y=mean),colour="gray25",linetype=3,size=0.5)+
  scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000),labels=c("0",expression(10^1),expression(10^2), expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7),expression(10^8),expression(10^9)))+
  #scale_y_continuous(expand= c(0,0)) +
  scale_x_continuous(expand= c(0,0)) +
  scale_color_brewer(palette='Set1') +
  scale_shape_manual(values=c(0,1,5))+
  ylab(ylabel)+
  xlab(xlabel)+
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.25, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0)),
          text = element_text(size=12),
          axis.title=element_text(size=14,face="bold"))

  return(g)
}