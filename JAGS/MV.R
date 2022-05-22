N_err<-GCS$N_GC_err
err_MV_T<-GCS$M_V_err
N = nrow(GCS)

######## NB with errors ########################################################
MV_Tx = seq(from = 1.05 * min(GCS$M_V), 
            to = 0.95 * max(GCS$M_V), 
            length.out = 500)

jags.data <- list(
  N_GC = GCS$N_GC,
  MV_T = GCS$M_V,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  err_MV_T = err_MV_T,
  MV_Tx = MV_Tx,
  M = 500
 )


model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)

# Prior for size

size~dunif(0.001,20)

#

for (i in 1:N){

MV_T_true[i]~dunif(-26,-10)
}

# Likelihood function

for (i in 1:N){

MV_T[i]~dnorm(MV_T_true[i],1/err_MV_T[i]^2);

errorN[i]~dbin(0.5,2*errN_GC[i])

eta[i]<-beta.0+beta.1*MV_T_true[i]

#log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems.
log(mu[i])<-log(exp(eta[i])+errorN[i]-errN_GC[i])

# Using dnegbin
N_GC[i]~dnegbin(p[i],size)
p[i]<-size/(size+mu[i])

#Using mixture of Poisson and Gamma
#N_GC[i]~dpois(g[i])
#g[i]~dgamma(size,rateParm[i])
#rateParm[i]<-size/mu[i]

# Prediction
etaTrue[i]<-beta.0+beta.1*MV_T_true[i]
    log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
    pTrue[i]<-size/(size+muTrue[i])
prediction.NB[i]~dnegbin(pTrue[i],size)
#prediction.NB[i]~dnegbin(p[i],size)

# Discrepancy measures
YNew[i] ~ dnegbin(p[i],size)
#YNew[i] ~ dpois(g[i])
expY[i] <- mu[i]
varY[i] <- mu[i] + pow(mu[i],2) / size
PRes[i] <-(N_GC[i] - expY[i])/sqrt(varY[i])
PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)

}
Fit<-sum(D[1:N])
New<-sum(DNew[1:N])

# Prediction for new data
for (j in 1:M){
  etax[j]<-beta.0+beta.1*MV_Tx[j]
  log(mux[j])<-max(-20,min(20,etax[j]))
  px[j]<-size/(size+mux[j])
  prediction.NBx[j]~dnegbin(px[j],size)
#prediction.NBx[j]~dpois(gx[j])
#gx[j]~dgamma(size,rateParmx[j])
#rateParmx[j]<-size/mux[j]
}
}"
inits1 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits2 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits3 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
params <- c("beta.0","beta.1","size","PRes","MV_T_true","Fit","New","prediction.NBx")

#jags.neg <- jags.model(
#  data = jags.data, 
#  inits = inits, 
#  textConnection(model.NB),
#  n.chains = 3,
#  n.adapt=1000
#)


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
summary<-extend.jags(jags.neg,drop.monitor=c("PRes","MV_T_true","Fit","New","prediction.NBx"), summarise=TRUE)

# Check the output
print(summary)

#update(jags.neg, 10000)

#jagssamples.nb <- jags.samples(jags.neg, params, n.iter = 50000)


#MV_T_true<-summary(as.mcmc.list(jags.neg,vars="MV_T_true"),quantiles=0.5)
#pred.NBerr<-summary(as.mcmc.list(jagssamples.nb, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
#pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MV_T_true=MV_T_true$quantiles,MV_T=GCS$MV_T,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])
pred.NBerrx<-summary(as.mcmc.list(jags.neg, vars="prediction.NBx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2errx<-data.frame(MV_Tx=MV_Tx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])



# N_low<-GCS$N_GC-N_err

# N_low[N_low<0]<-0


# GCS$Type <- factor(GCS$group_type, levels = c("E", "Spiral", "S0", "Irr"))
# g2<-ggplot(GCS,aes(x=M_V,y=N_GC))+
#   geom_ribbon(data=pred.NB2errx,aes(x=MV_Tx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray") +
#   geom_ribbon(data=pred.NB2errx,aes(x=MV_Tx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill="gray") +
#   geom_ribbon(data=pred.NB2errx,aes(x=MV_Tx,y=mean,ymin=lwr3, ymax=upr3), alpha=0.25, fill="gray") +
#   geom_point(aes(colour=Type,shape=Type),size=3.25,alpha=0.8)+
#   geom_errorbar(aes(colour=Type,ymin=N_low,ymax=N_GC+N_err),alpha=0.7,width=0.05)+
#   geom_errorbarh(aes(colour=Type,xmin=M_V-GCS$M_V_err, xmax=M_V+M_V_err),alpha=0.7,height=0.05)+
#   geom_line(data=pred.NB2errx,aes(x=MV_Tx,y=mean),colour="gray25",linetype="dashed",size=1.2)+
#   scale_y_continuous(trans='asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
#   scale_colour_gdocs()+
#   scale_shape_manual(values=c(19,2,8,10))+scale_x_reverse()+
#   theme_hc()+
#   ylab(ylabel)+
#   xlab(expression(M[V]))+theme(legend.position="top",plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25),text = element_text(size=25))

# CairoPDF(paste(datapath,'M_V.pdf',sep=''),height=8,width=9)
# g2
# dev.off()


# # Model comparison 

# Pred<-ggs(jagssamples.nb,family=c("New"))[,"value"]
# Obs<-ggs(jagssamples.nb,family=c("Fit"))[,"value"]
# sqrt(mean((Pred-Obs)^2))

# # Dispersion parameter

# require(scales)
Pres<-summary(as.mcmc.list(jags.neg, vars="PRes"),quantiles=0.5)$quantiles
Dispersion = sum(Pres^2)/(N-3)# beta.0, beta.1 and k, 3 parameters

# DIC 

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
