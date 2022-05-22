N = nrow(df)

######## Gamma GLM with errors in both axis ########################################################
dfx = seq(from = 0.95 * min(df$x), 
           to = 1.05 * max(df$x), 
           length.out = 500)

dfz = seq(from = 0.95 * min(df$z), 
           to = 1.05 * max(df$z), 
           length.out = 500)

## Define dataset for JAGS#######################################

jags.data <- list(
  y = df$y,
  x = df$x,
  z = df$z,
  erry = df$y_err,
  N = nrow(df),
  errx = df$up_x,
  errz = df$up_z,
  dfx = dfx,
  dfz = dfz,
  M = 500
)

#### JAGS model#######################################
model.gamma <- "model{

### Priors

# Regression coefficients - noninformative

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
beta.2~dnorm(0,0.000001)

#shape - noninformative
shape ~ dunif(-100,100)

#Hyperpriors
meanx ~ dgamma(0.01,0.01)
varx ~ dgamma(0.01,0.01)

meanz ~ dgamma(0.01,0.01)
varz ~ dgamma(0.01,0.01)

# 1.Likelihood function

for (i in 1:N){

y[i]~dgamma(shape, shape/mu[i])
eta[i]<-beta.0+beta.1*x_true[i]+beta.2*z_true[i]
log(mu[i])<-log(exp(eta[i])+errorY[i])
errorY[i]~dnorm(0, 1/(erry[i]^2))

x_true[i] ~ dgamma(meanx^2/varx,meanx/varx)# True predictor
x[i]~dnorm(x_true[i],1/(errx[i]^2))# Observed predictor
z_true[i] ~ dgamma(meanz^2/varz,meanz/varz)# True predictor
z[i]~dnorm(z_true[i],1/(errz[i]^2))# Observed predictor

# # 2.Discrepancy measures

# YNew[i] ~ dnorm(mu[i],tau)
# #expY[i] <- mu[i]
# varY[i] <- 1 / tau
# PRes[i] <-(y[i] - mu[i])/sqrt(varY[i])
# PResNew[i] <-(YNew[i] - mu[i])/sqrt(varY[i])
# D[i]<-pow(PRes[i],2)
# DNew[i]<-pow(PResNew[i],2)

# # 3.Prediction of current data

  # etaTrue[i]<-beta.0+beta.1*x[i]
  # log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
  # prediction.gamma[i]~dgamma(muTrue[i],shape)
}
# Fit<-sum(D[1:N])
# New<-sum(DNew[1:N])

# # 4.Prediction new data

# for (j in 1:M){
#   mux[j]<-beta.0+beta.1*dfx[j]
# #  log(mux[j])<-max(-20,min(20,etax[j]))
#   prediction.Normx[j]~dnorm(mux[j],tau)
# }
}"

# inits1 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
# inits2 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
# inits3 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
# params <- c("beta.0","beta.1","tau","PRes","x_true","Fit","New","prediction.Normx")

inits<-function(){list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),beta.2=rnorm(1,0,0.1),shape=runif(1,0.1,5))}

# One for each chain
inits1 <- inits()
inits2 <- inits()
inits3 <- inits()
#params <- c("beta.0","beta.1","size","PRes","prediction.NB","xtrue","Fit","New","prediction.NBx")
params <- c("beta.0","beta.1","beta.2","shape",'dic')

#monparams <- c("beta.0","beta.1","beta.2","size

#jags.neg <- jags.model(
#  data = jags.data, 
#  inits = inits, 
#  textConnection(model.NB),
#  n.chains = 3,
#  n.adapt=1000
#)

jags.gamma <- autorun.jags(method="rjags",
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model.gamma,
                     n.chains = 3,
                     adapt=2000,
                     monitor=c(params),
                     startburnin=burnin,
                     startsample=sample,
                     thin=thin,
                     summarise=FALSE,
                     plots=FALSE
)

summary<-extend.jags(jags.gamma, summarise=TRUE)
#dicsample<-extend.jags(jags.gamma, drop.monitor=params, add.monitor=c('dic'))
#jagssamples.g3 <- as.mcmc.list(summary)

# # Check the output
print(summary)

# #update(jags.neg, 10000)

# jags.gamma <- jags.model(
#       data = jags.data,
#       inits=inits3,
#       textConnection(model.gamma),
#       n.chains=3,
#       n.adapt=5000
# )

# #jagssamples.gamma <- jags.samples(jags.gamma, n.iter = 50000)
# update(jags.gamma, 100000)

# jagssamples.g3 <- jags.samples(jags.gamma, params, n.iter = 5000)
# codasamples.g3 <- coda.samples(jags.gamma, params, n.iter = 5000)
# dicsamples.g3 <- dic.samples(jags.gamma, params, n.iter = 5000,type="pD")

# #mass_true<-summary(as.mcmc.list(jags.neg,vars="mass_true"),quantiles=0.5)
#pred.Gerr<-summary(as.mcmc.list(jagssamples.g3, vars="prediction.gamma"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# #pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MV_T_true=MV_T_true$quantiles,MV_T=GCS$MV_T,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])
# pred.Normerrx<-summary(as.mcmc.list(jags.norm, vars="prediction.Normx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# pred.Norm2errx<-data.frame(x=x2,mean=pred.Normerrx$statistics[,1],lwr1=pred.Normerrx$quantiles[,3],lwr2=pred.Normerrx$quantiles[,2],lwr3=pred.Normerrx$quantiles[,1],upr1=pred.Normerrx$quantiles[,5],upr2=pred.Normerrx$quantiles[,6],upr3=pred.Normerrx$quantiles[,7])

# S.G1<-ggs(jagssamples.g3,family=c("beta"))
# S.G2<-ggs(jagssamples.g3,family=c("shape"))

# S.G<-rbind(S.G1,S.G2,deparse.level=2)
# ggs_density(S.G)+theme_fivethirtyeight()

# N_low<-GCS$Extent-N_err

# N_low[N_low<0]<-0


# GCS$Type <- factor(GCS$group_type, levels = c("E", "Spiral", "S0", "Irr"))
# g2<-ggplot(GCS,aes(x=Mdyn,y=Extent))+
#   geom_ribbon(data=pred.Norm2errx,aes(x=x2,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray") +
#   geom_ribbon(data=pred.Norm2errx,aes(x=x2,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill="gray") +
#   geom_ribbon(data=pred.Norm2errx,aes(x=x2,y=mean,ymin=lwr3, ymax=upr3), alpha=0.25, fill="gray") +
#   geom_point(aes(colour=Type,shape=Type),size=3.25,alpha=0.8)+
#   geom_errorbar(aes(colour=Type,ymin=N_low,ymax=Extent+N_err),alpha=0.7,width=0.05)+
#   geom_errorbarh(aes(colour=Type,xmin=Mdyn-Mdyn_err, xmax=Mdyn+Mdyn_err),alpha=0.7,height=0.05)+
#   geom_line(data=pred.Norm2errx,aes(x=x2,y=mean),colour="gray25",linetype="dashed",size=0.8)+
#   scale_y_continuous(trans='asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
#   scale_colour_gdocs()+
#   scale_shape_manual(values=c(19,2,8,10))+
#   theme_hc()+
#   theme(legend.position="top", plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25), text = element_text(size=25)) +
#   ylab("Extent")+
#   xlab(expression(log~M[dyn]/M['\u0298']))

# CairoPDF('Mdyn.pdf',height=8,width=9)
# g2
# dev.off()



# Pred<-ggs(jagssamples.nb,family=c("New"))[,"value"]
# Obs<-ggs(jagssamples.nb,family=c("Fit"))[,"value"]
# sqrt(mean((Pred-Obs)^2))

# # Dispersion parameter

# # require(scales)
# Pres<-summary(as.mcmc.list(jags.norm, vars="PRes"),quantiles=0.5)$quantiles
# Dispersion = sum(Pres^2)/(N-3)# beta.0, beta.1 and k, 3 parameters

dicinits <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
jags.DIC <- jags.model(
  data = jags.data, 
  inits = inits1, 
  textConnection(model.gamma),
  n.chains = 3,
  n.adapt=2000
)

# update(jags.DIC , 10000)
# dicsamples.norm <- dic.samples(jags.DIC, params, n.iter = 25000,type="pD")

# dicsamples.norm











