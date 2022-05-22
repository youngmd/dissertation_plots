N = nrow(df)

######## NB GLM with errors in both axis ########################################################

### Vector of new values for prediction
dfx = seq(from = 0.95 * min(df$x), 
           to = 1.05 * max(df$x), 
           length.out = 500)


## Define dataset for JAGS#######################################

jags.data <- list(
  y = df$y,
  x = df$x,
  erry = df$y_err,
  N = nrow(df),
  errx = df$up_x,
  dfx = dfx,
  M = 500
)

#### JAGS model#######################################
model.lnorm <- "model{

### Priors

# Regression coefficients - noninformative

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)

#tau - noninformative
tau ~ dgamma(0.01, 0.01)
sigma <- 1/sqrt(tau)

#Hyperpriors

meanx ~ dgamma(0.01,0.01)
varx ~ dgamma(0.01,0.01)

# 1.Likelihood function

for (i in 1:N){

y[i]~dlnorm(mu[i],tau)
mu[i]<- beta.0+beta.1*x_true[i] + errorY[i]

errorY[i]~dnorm(0, 1/(erry[i]^2))

x_true[i] ~ dgamma(meanx^2/varx,meanx/varx)# True predictor
x[i]~dnorm(x_true[i],1/(errx[i]^2))# Observed predictor


# # 2.Discrepancy measures

# YNew[i] ~ dnorm(mu[i],tau)
# #expY[i] <- mu[i]
# varY[i] <- 1 / tau
# PRes[i] <-(y[i] - mu[i])/sqrt(varY[i])
# PResNew[i] <-(YNew[i] - mu[i])/sqrt(varY[i])
# D[i]<-pow(PRes[i],2)
# DNew[i]<-pow(PResNew[i],2)

# # 3.Prediction of current data

  # muTrue[i]<-exp(beta.0+beta.1*x_true[i])
  # # log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
  # prediction.Norm[i]~dnorm(muTrue[i],tau)
}
# Fit<-sum(D[1:N])
# New<-sum(DNew[1:N])

# # 4.Prediction new data

for (j in 1:M){
  mux[j]<-exp(beta.0+beta.1*dfx[j])
  # log(mux[j])<-max(-20,min(20,etax[j]))
  prediction.lNormx[j]~dlnorm(mux[j],tau)
}
}"

# inits1 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
# inits2 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
# inits3 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),sigma=runif(1,0.1,5))
# params <- c("beta.0","beta.1","tau","PRes","x_true","Fit","New","prediction.Normx")

inits<-function(){list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),tau=runif(1,0.1,20))}

# One for each chain
inits1 <- inits()
inits2 <- inits()
inits3 <- inits()
#params <- c("beta.0","beta.1","size","PRes","prediction.NB","xtrue","Fit","New","prediction.NBx")
params <- c("beta.0","beta.1","sigma")

#monparams <- c("beta.0","beta.1","beta.2","size

#jags.neg <- jags.model(
#  data = jags.data, 
#  inits = inits, 
#  textConnection(model.NB),
#  n.chains = 3,
#  n.adapt=1000
#)

jags.lnorm <- autorun.jags(method="rjparallel",
                     method.options=list(cl=cl),
                     modules=c('glm'),
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model.lnorm,
                     n.chains = 3,
                     adapt=2000,
                     monitor=c(params),
                     startburnin=burnin,
                     startsample=sample,
                     thin=thin,
                     summarise=FALSE,
                     plots=FALSE,
)

summary<-extend.jags(jags.lnorm, summarise=TRUE)
ssamples.n <- as.mcmc.list(summary)

S.G1<-ggs(ssamples.n,family=c("beta"))
S.G2<-ggs(ssamples.n,family=c("sigma"))

S.G<-rbind(S.G1,S.G2,deparse.level=2)

postplots <- ggs_density(S.G)+theme_fivethirtyeight()
traceplots <- ggs_traceplot(S.G)+theme_fivethirtyeight()

jagssamples.lnorm<-extend.jags(jags.norm, drop.monitor=c('beta.0','beta.1','sigma'), add.monitor=c('prediction.lNormx'), summarise=FALSE)
# summary<-extend.jags(jags.norm,drop.monitor=c("PRes","x_true","Fit","New","prediction.Normx"), summarise=TRUE)

# # Check the output
print(summary)

# #update(jags.neg, 10000)

# #jagssamples.nb <- jags.samples(jags.neg, params, n.iter = 50000)


# #mass_true<-summary(as.mcmc.list(jags.neg,vars="mass_true"),quantiles=0.5)
# #pred.NBerr<-summary(as.mcmc.list(jagssamples.nb, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# #pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MV_T_true=MV_T_true$quantiles,MV_T=GCS$MV_T,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])
pred.Normerrx<-summary(as.mcmc.list(jagssamples.lnorm, vars="prediction.lNormx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.Norm2errx<-data.frame(x=dfx,mean=pred.Normerrx$statistics[,1],lwr1=pred.Normerrx$quantiles[,3],lwr2=pred.Normerrx$quantiles[,2],lwr3=pred.Normerrx$quantiles[,1],upr1=pred.Normerrx$quantiles[,5],upr2=pred.Normerrx$quantiles[,6],upr3=pred.Normerrx$quantiles[,7])



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

# CairoPDF('norm.pdf',height=8,width=9)
# g2
# dev.off()



# Pred<-ggs(jagssamples.nb,family=c("New"))[,"value"]
# Obs<-ggs(jagssamples.nb,family=c("Fit"))[,"value"]
# sqrt(mean((Pred-Obs)^2))

# Dispersion parameter

# require(scales)
# Pres<-summary(as.mcmc.list(jags.norm, vars="PRes"),quantiles=0.5)$quantiles
# Dispersion = sum(Pres^2)/(N-3)# beta.0, beta.1 and k, 3 parameters

# #use the mean of the end results to initiate the DIC sample
inits <- list('beta.0'=summary$summaries['beta.0','Mean'], 'beta.1'=summary$summaries['beta.1','Mean'],tau=1/(sqrt(summary$summaries['sigma','Mean'])))

jags.DIC <- jags.model(
  data = jags.data, 
  inits = inits, 
  textConnection(model.lnorm),
  n.chains = 3,
  n.adapt=2000
)

update(jags.DIC, 30000)
dicsamples.lnorm <- dic.samples(jags.DIC, params, n.iter = 15000,type="pD")

dicsamples.lnorm











