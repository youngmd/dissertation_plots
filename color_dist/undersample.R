library(mixtools)
library(mclust)
library(ggplot2)
library(Cairo)

gcc90<- read.csv('gcc90.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

gccdata <- subset(gcc90, galaxy_name == 'NGC5846')

makesample <- function(N, foo) {
  z <- rbinom(N, 2, foo$lambda)
  z1 <- z
  z1[which(z==0)] = 1
  z1[which(z!=0)] = 0
  z2 <- z
  z2[which(z!=1)] = 0
  z3 <- z
  z3[which(z==2)] = 1
  z3[which(z!=2)] = 0
  x.b <- z1*rnorm(N, foo$mu[1], foo$sigma[1]) + z2*rnorm(N, foo$mu[2], foo$sigma[2]) + z3*rnorm(N, foo$mu[3], foo$sigma[3])
  return(x.b)
}

n <- nrow(gccdata)
br <- gccdata$BR

df <- data.frame('n'=0,'avg'=0, 'std'=0)

iters = 500
m <- Mclust(br, G=3, modelName='E')

m3 <- normalmixEM(br, mu=m$parameters$mean, maxit=1000, lambda=m$parameters$pro, sigma=sqrt(m$parameters$variance$sigmasq), k=3, mean.constr=c(NA, NA, NA), arbmean=TRUE, arbvar=FALSE)

n = 1500

for( i in seq(20,n,20)){
    print(i)
    tmp <- rep(0,iters)
    for( j in 1:iters){
        ss <- makesample(n-i, m3)
        m <- Mclust(ss, modelName="E")
        best <- m$bic
        ebics <- m$BIC[,'E']
        ebics <- subset(ebics, abs(ebics - best) < 5)
        G <- mean(as.numeric(names(ebics)))
        tmp[j] <- G
    }
    av <- mean(tmp)
    st <- sd(tmp)
    dft <- data.frame('n'=n-i,'avg'=av,'std'=st)
    df <- rbind(df, dft)
    if(n-i < 30){
        break
    }
}

df <- subset(df, n != 0)
df <- subset(df, n > 30)

d2 <- read.table('mm_comp.dat', header=TRUE)

d2$Type <- factor(d2$type, levels = c("E", "Spiral", "S0"))
g <- ggplot(data=df, aes(x=n, y=avg)) + 
  geom_ribbon(aes(ymin=avg-std, ymax=avg+std), fill="grey40", alpha=0.2) + 
  geom_line(linetype=2, color="black", size=0.5) + 
  geom_point(data=d2, aes(x=n90, y=nmm, shape=Type, color=Type), size=3) + 
  scale_x_log10(expand = c(0,0)) +
  theme_bw() + 
  theme(legend.position=c(0.7,0.2)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        legend.key = element_blank(),
        axis.text.y=element_text(margin=margin(0,15,0,0)), 
        axis.text.x=element_text(margin=margin(15,0,0,0))) +
  xlab(expression('Sample Size of GC Candidate B-R Color Distribution')) + 
  scale_shape_manual(values=c(2,3,4))+
  ylab(expression('Gaussian Mixture Model Components')) +
  scale_y_continuous(limits=c(0.5,3.5)) + 
  annotation_logticks(sides = "b")

ggsave("mm_trend.pdf", g, width=8, height=8)