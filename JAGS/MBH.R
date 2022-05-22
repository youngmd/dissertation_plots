N_err<-GCS$N_GC_err
GCS$lowMBH<-GCS$lg_M_B_err_low
GCS$upMBH<-GCS$lg_M_B_err_high
N = nrow(GCS)

######## NB GLM with errors in both axis ########################################################

### Vector of new values for prediction
MBHx = seq(from = 0.95 * min(GCS$lg_M_B), 
           to = 1.05 * max(GCS$lg_M_B), 
           length.out = 500)


## Define dataset for JAGS#######################################

jags.data <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$lg_M_B,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  errMBH = GCS$upMBH,
  MBHx = MBHx,
  M = 500
)


#### JAGS model#######################################
model.NB <- "model{

### Priors

# Regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)

#Size

size~dunif(0.001,20)

#Hyperpriors

meanx ~ dgamma(0.01,0.01)
varx ~ dgamma(0.01,0.01)

# 1.Likelihood function

for (i in 1:N){

N_GC[i]~dnegbin(p[i],size)

p[i]<-size/(size+mu[i])

eta[i]<-beta.0+beta.1*MBHtrue[i]

log(mu[i])<-log(exp(eta[i])+errorN[i]-errN_GC[i])

errorN[i]~dbin(0.5,2*errN_GC[i])# Errors in N_GC

MBHtrue[i] ~ dgamma(meanx^2/varx,meanx/varx)# True black hole mass (predictor)

MBH[i]~dnorm(MBHtrue[i],1/errMBH[i]^2)# Observed black hole mass


# 2.Discrepancy measures

YNew[i] ~ dnegbin(p[i],size)
expY[i] <- mu[i]
varY[i] <- mu[i] + pow(mu[i],2) / size
PRes[i] <-(N_GC[i] - expY[i])/sqrt(varY[i])
PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)

# 3.Prediction of current data

  etaTrue[i]<-beta.0+beta.1*MBHtrue[i]
  log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
  pTrue[i]<-size/(size+muTrue[i])
  prediction.NB[i]~dnegbin(pTrue[i],size)
}
Fit<-sum(D[1:N])
New<-sum(DNew[1:N])

# 4.Prediction new data

for (j in 1:M){
  etax[j]<-beta.0+beta.1*MBHx[j]
  log(mux[j])<-max(-20,min(20,etax[j]))
  px[j]<-size/(size+mux[j])
  prediction.NBx[j]~dnegbin(px[j],size)
}
}"
############  

# Define function to generate initial values
inits<-function(){list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))}

# One for each chain
inits1 <- inits()
inits2 <- inits()
inits3 <- inits()
params <- c("beta.0","beta.1","size","PRes","prediction.NB","MBHtrue","Fit","New","prediction.NBx")



# Run jags in parallel 
jags.neg <- run.jags(method="rjparallel", method.options=list(cl=cl),
  data = jags.data, 
  inits = list(inits1,inits2,inits3),
  model=model.NB,
  n.chains = 3,
  adapt=2000,
  monitor=c(params),
   burnin=burnin,
  sample=sample,
  thin=thin,
  summarise=FALSE,
  plots=FALSE
)

jagssamples.nb <- as.mcmc.list(jags.neg )
summary<-extend.jags(jags.neg,drop.monitor=c("PRes","prediction.NB","MBHtrue","Fit","New","prediction.NBx"), summarise=TRUE)

# Check the output
print(summary)

## Analysis ends here####


print('starting diagnostics')
# ## Below are plot diagnostics
# MBHtrue<-summary(as.mcmc.list(jags.neg, vars="MBHtrue"),quantiles=0.5)
# pred.NBerr<-summary(as.mcmc.list(jags.neg, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# pred.NB2err<-data.frame(Type=GCS$type,NGC=GCS$Nblue,MBHtrue=MBHtrue$quantiles,MBH=GCS$M_BH,mean=pred.NBerr$statistics[1:23,1],lwr1=pred.NBerr$quantiles[1:23,3],lwr2=pred.NBerr$quantiles[1:23,2],lwr3=pred.NBerr$quantiles[1:23,1],upr1=pred.NBerr$quantiles[1:23,5],upr2=pred.NBerr$quantiles[1:23,6],upr3=pred.NBerr$quantiles[1:23,7])
pred.NBerrx<-summary(as.mcmc.list(jags.neg,vars="prediction.NBx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2errx<-data.frame(MBHx=MBHx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])

print('done with diagnostics')

# #N_low<-asinh(pred.NB2err$NGC-N_err)
# N_low<-GCS$N_GC-N_err
# N_low[N_low<0]<-0

# asinh_trans <- function(){
#   trans_new(name = 'asinh', transform = function(x) asinh(x), 
#             inverse = function(x) sinh(x))
# }

# print('Starting plot')
# GCS$Type <- factor(GCS$group_type, levels = c("E", "Spiral", "S0"))
# g2 <- ggplot(GCS,aes(x=lg_M_B,y=N_GC))+
#   geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.35, fill="gray") +
#   geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.25, fill="gray") +
#   geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr3, ymax=upr3), alpha=0.15, fill="gray") +
#   geom_point(aes(x=lg_M_B, y=N_GC, colour=Type,shape=Type),size=2.5,alpha=0.8)+
#   geom_errorbar(aes(colour=Type,ymin=N_low,ymax=N_GC+N_err),alpha=0.7,width=0.03)+
#   geom_errorbarh(aes(colour=Type,xmin=lg_M_B-lowMBH, xmax=lg_M_B+upMBH),alpha=0.7,height=0.03)+
#   geom_line(data=pred.NB2errx,aes(x=MBHx,y=mean),colour="gray25",linetype="dashed",size=0.6)+
#   scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2), expression(10^3),expression(10^4),expression(10^5)))+
#   scale_colour_economist()+
#   scale_shape_manual(values=c(2,3,4))+
#   theme_hc()+
#   ylab(ylabel)+
#   xlab(expression(log~M[BH]/M['\u0298']))+
#   theme(legend.position="none",plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25), text = element_text(size=20), panel.grid.major = element_blank())

# CairoPDF(paste(datapath,'mass_bh.pdf',sep=''),height=8,width=9)
# g2
# dev.off()


# Pred<-ggs(jagssamples.nb,family=c("New"))[,"value"]
# Obs<-ggs(jagssamples.nb,family=c("Fit"))[,"value"]
# sqrt(mean((Pred-Obs)^2))

# Dispersion parameter

# require(scales)
Pres<-summary(as.mcmc.list(jags.neg, vars="PRes"),quantiles=0.5)$quantiles
Dispersion = sum(Pres^2)/(N-3)# beta.0, beta.1 and k, 3 parameters


jags.DIC <- jags.model(
  data = jags.data, 
  inits = inits1, 
  textConnection(model.NB),
  n.chains = 3,
  n.adapt=2000
)

update(jags.DIC , 10000)
dicsamples.nb <- dic.samples(jags.DIC, params, n.iter = 25000,type="pD")

dicsamples.nb





