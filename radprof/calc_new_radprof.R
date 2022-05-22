# Bayesian NB regression using JAGS by 
#Rafael S. de Souza, Bart Buelens, Ewan Cameron, Joseph Hilbe

#########################################
##Required libraries#####################
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
library(runjags)
library(parallel)

##########################################


###########################################
##Auxiliar functions for ploting 
# Function to allow parse labels in facet_wrap

facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)
  
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}
give.n <- function(x){
  
  return(c(y = 0.5, label = length(x))) 
  #
}

# Arcsinh transformation
asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

##################################################################################


# Define working directory if necessary
#dir<-getwd()
#setwd(dir)


###Script starts here#################

##Read data#################
cnames <- c("r","N","sd","sderr","lsd","lsderrhi","lsderrlo","good_frac","r_d","NGC","Ncontam","contam_frac")

rp = read.table(paste(datapath,'/corrected_profile.all.out',sep=''), col.names=cnames)

#GCS = read.csv(file="GCs.csv",header=TRUE,dec=".",sep="")
rp = subset(rp, sd > 0) # 1 removed
sd_err<-rp$sderr
N = nrow(rp)

######## NB GLM with errors in both axis ########################################################

### Vector of new values for prediction
Rx = seq(from = 0.0001, 
           to = 1.05 * max(rp$r), 
           length.out = 500)


## Define dataset for JAGS#######################################

jags.data <- list(
  R = rp$r,
  SD = rp$sd,
  errSD = rp$sderr,
  N = nrow(rp),
  Rx = Rx,
  M = 500
)


#### JAGS model#######################################
model.Pois <- "model{

### Priors

# Regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)

# 1.Likelihood function

for (i in 1:N){

SD[i]~dpois(mu[i])

eta[i]<-beta.0+beta.1*log(R[i])

log(mu[i])<- eta[i]


# 2.Discrepancy measures

YNew[i] ~ dpois(mu[i])
expY[i] <- mu[i]
varY[i] <- mu[i] + pow(mu[i],2)
PRes[i] <-(SD[i] - expY[i])/sqrt(varY[i])
PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)

# 3.Prediction of current data

  etaTrue[i]<-beta.0+beta.1*log(R[i])
  log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
  prediction.P[i]~dpois(muTrue[i])
}
Fit<-sum(D[1:N])
New<-sum(DNew[1:N])

# 4.Prediction new data

for (j in 1:M){
  etax[j]<-beta.0+beta.1*log(Rx[j])
  log(mux[j])<-max(-20,min(20,etax[j]))
  prediction.Px[j]~dpois(mux[j])
}
}"
############  

# Define function to generate initial values
inits<-function(){list(beta.0=rnorm(1,3,1),beta.1=rnorm(1,-3,1))}

# One for each chain
inits1 <- inits()
inits2 <- inits()
inits3 <- inits()
params <- c("beta.0","beta.1","PRes","prediction.P","Fit","New","prediction.Px")



# Allocate 3 cores in the computer
cl <- makeCluster(8)
##

# Run jags in parallel 
jags.neg <- run.jags(method="rjparallel", method.options=list(cl=cl),
  data = jags.data, 
  inits = list(inits1,inits2,inits3),
  model=model.Pois,
  n.chains = 3,
  adapt=2000,
  monitor=c(params),
   burnin=20000,
  sample=50000,
  summarise=FALSE,
  plots=FALSE
)

jagssamples.nb <- as.mcmc.list(jags.neg )
summary<-extend.jags(jags.neg,drop.monitor=c("PRes","prediction.P","Fit","New","prediction.Px"), summarise=TRUE)

# Check the output
print(summary)

### Analysis ends here####


# print('starting diagnostics')
# ## Below are plot diagnostics
# MBHtrue<-summary(as.mcmc.list(jags.neg, vars="MBHtrue"),quantiles=0.5)
# pred.NBerr<-summary(as.mcmc.list(jags.neg, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MBHtrue=MBHtrue$quantiles,MBH=GCS$MBH,mean=pred.NBerr$statistics[1:45,1],lwr1=pred.NBerr$quantiles[1:45,3],lwr2=pred.NBerr$quantiles[1:45,2],lwr3=pred.NBerr$quantiles[1:45,1],upr1=pred.NBerr$quantiles[1:45,5],upr2=pred.NBerr$quantiles[1:45,6],upr3=pred.NBerr$quantiles[1:45,7])
# pred.NBerrx<-summary(as.mcmc.list(jags.neg,vars="prediction.NBx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# pred.NB2errx<-data.frame(MBHx=MBHx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])

# print('done with diagnostics')

# #N_low<-asinh(pred.NB2err$NGC-N_err)
# N_low<-pred.NB2err$NGC-N_err
# N_low[N_low<0]<-0

# print('Starting plot')
# GCS$Type <- factor(GCS$Type, levels = c("E", "S", "S0", "Irr"))
# cairo_pdf("MBHx2.pdf",height=8,width=9)
# ggplot(pred.NB2err,aes(x=MBH,y=NGC))+
#   geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray") +
#   geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill="gray") +
#   geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr3, ymax=upr3), alpha=0.25, fill="gray") +
#   geom_point(aes(colour=Type,shape=Type),size=3.25,alpha=0.8)+
#   geom_errorbar(guide="none",aes(colour=Type,ymin=N_low,ymax=NGC+N_err),alpha=0.7,width=0.05)+
#   geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
#                                   xmax=MBH+upMBH),alpha=0.7,height=0.05)+
#   geom_line(data=pred.NB2errx,aes(x=MBHx,y=mean),colour="gray25",linetype="dashed",size=1.2)+
#   annotate("text", x = 6.63, y = 800, label = "Milky Way",size = 6.5)+
#   geom_segment(aes(x =  6.65, y = 600, xend = 6.61, yend = 200), arrow = arrow(length = unit(0.25, "cm")))+
#   scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),
#                                                                                    expression(10^3),expression(10^4),expression(10^5)))+
#   scale_colour_gdocs()+
#   scale_shape_manual(values=c(19,2,8))+
#   theme_hc()+
#   ylab(expression(N[GC]))+
#   xlab(expression(log~M[BH]/M['\u0298']))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
#                                                 axis.title.y=element_text(vjust=0.75),
#                                                 axis.title.x=element_text(vjust=-0.25),
#                                                 text = element_text(size=25))
# dev.off()

ggplot(data=rp, (aes(r, sd))+geom_point()+geom_errorbar(aes(ymin=sd-sderr,ymax=sd+sderr), alpha=0.7,width=0.05)+geom_smooth(method='glm', aes(x=log(r), y=log(sd)), level=0.9)