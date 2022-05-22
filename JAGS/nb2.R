# y_err<-df$y_err
# df$downx<-df$xlow
# df$upx<-df$xhigh
N = nrow(df)

######## NB GLM with errors in both axis ########################################################

### Vector of new values for prediction
# dfx = seq(from = 0.95 * min(df$x), 
#            to = 1.05 * max(df$x), 
#            length.out = 500)
# 
# dfz = seq(from = 0.95 * min(df$z), 
#            to = 1.05 * max(df$z), 
#            length.out = 500)

## Define dataset for JAGS#######################################

jags.data <- list(
  y = df$y,
  x = df$x,
  z = df$z,
  erry = df$y_err,
  N = nrow(df),
  errx = df$up_x,
  errz = df$up_z,
  xtmin = 0.8*min(df$x),
  xtmax = 1.2*max(df$x),
  ztmin = 0.8*min(df$z),
  ztmax = 1.2*max(df$z)
)


#### JAGS model#######################################
model.NB2 <- "model{

### Priors

# Regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
beta.2~dnorm(0,0.000001)

#Size

size~dunif(0.001,50)

# #Hyperpriors
# 
# meanx ~ dgamma(0.01,0.01)
# varx ~ dgamma(0.01,0.01)
# 
# meanz ~ dgamma(0.01,0.01)
# varz ~ dgamma(0.01,0.01)

for (i in 1:N){
  xtrue[i] ~ dunif(xtmin,xtmax)# True predictor1
  ztrue[i] ~ dunif(ztmin,ztmax)# True predictor2
}

# 1.Likelihood function

for (i in 1:N){

x[i]~dnorm(xtrue[i],1/errx[i]^2)# Observed predictor1
z[i]~dnorm(ztrue[i],1/errz[i]^2)# Observed predictor2

y[i]~dnegbin(p[i],size)
p[i]<-size/(size+mu[i])
eta[i]<-beta.0+beta.1*xtrue[i] + beta.2*ztrue[i] + exp(errorY[i]-erry[i])
log(mu[i])<-max(-20,min(20,eta[i]))
errorY[i]~dbin(0.5,2*erry[i])# Errors in y


# # 2.Discrepancy measures

# YNew[i] ~ dnegbin(p[i],size)
# expY[i] <- mu[i]
# varY[i] <- mu[i] + pow(mu[i],2) / size
# PRes[i] <-(y[i] - expY[i])/sqrt(varY[i])
# PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
# D[i]<-pow(PRes[i],2)
# DNew[i]<-pow(PResNew[i],2)

# # 3.Prediction of current data

   etaTrue[i]<-beta.0+beta.1*xtrue[i]+beta.2*ztrue[i]
   log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
   pTrue[i]<-size/(size+muTrue[i])
   prediction.NB[i]~dnegbin(pTrue[i],size)
}
# Fit<-sum(D[1:N])
# New<-sum(DNew[1:N])

# # 4.Prediction new data

# for (j in 1:M){
#   etax[j]<-beta.0+beta.1*dfx[j]+beta.2*dfz[j]
#   log(mux[j])<-max(-20,min(20,etax[j]))
#   px[j]<-size/(size+mux[j])
#   prediction.NBx[j]~dnegbin(px[j],size)
# }
}"
############  

# Define function to generate initial values
inits<-function(){list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),beta.2=rnorm(1,0,0.1),size=runif(1,0.1,50))}

# One for each chain
inits1 <- inits()
inits2 <- inits()
inits3 <- inits()
params <- c("beta.0","beta.1","beta.2","size")

# Run jags in parallel 
jags.neg <- autorun.jags(method="rjparallel",
  method.options=list(cl=cl),
  modules=c('glm'),
  data = jags.data, 
  inits = list(inits1,inits2,inits3),
  model=model.NB2,
  n.chains = 3,
  adapt=5000,
  monitor=params,
  startburnin=burnin,
  startsample=sample,
  thin=thin,
  summarise=FALSE,
  plots=FALSE
)

summary<-extend.jags(jags.neg, summarise=TRUE)
summary

samples.nb <- extend.jags(jags.neg, drop.monitor=c("beta.0","beta.1","size"), add.monitor=c("prediction.NB"))

pred.NBerr<-summary(as.mcmc.list(samples.nb, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2err<-data.frame(Galaxy=df$galaxy,NGC=df$y,NGC_err=df$y_err,LK=df$x,Mdyn=df$z,mean=pred.NBerr$statistics[1:N,1],lwr1=pred.NBerr$quantiles[1:N,3],lwr2=pred.NBerr$quantiles[1:N,2],lwr3=pred.NBerr$quantiles[1:N,1],upr1=pred.NBerr$quantiles[1:N,5],upr2=pred.NBerr$quantiles[1:N,6],upr3=pred.NBerr$quantiles[1:N,7])


# print('starting diagnostics')
# # ## Below are plot diagnostics
# # MBHtrue<-summary(as.mcmc.list(jags.neg, vars="MBHtrue"),quantiles=0.5)
# # pred.NBerr<-summary(as.mcmc.list(jags.neg, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# # pred.NB2err<-data.frame(Type=df$type,NGC=df$Nblue,MBHtrue=MBHtrue$quantiles,MBH=df$M_BH,mean=pred.NBerr$statistics[1:23,1],lwr1=pred.NBerr$quantiles[1:23,3],lwr2=pred.NBerr$quantiles[1:23,2],lwr3=pred.NBerr$quantiles[1:23,1],upr1=pred.NBerr$quantiles[1:23,5],upr2=pred.NBerr$quantiles[1:23,6],upr3=pred.NBerr$quantiles[1:23,7])
# pred.NBerrx<-summary(as.mcmc.list(jags.neg,vars="prediction.NBx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
# pred.NB2errx<-data.frame(dfx=dfx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])

# print('done with diagnostics')

# #N_low<-asinh(pred.NB2err$NGC-N_err)
# N_low<-df$N_GC-N_err
# N_low[N_low<0]<-0

# asinh_trans <- function(){
#   trans_new(name = 'asinh', transform = function(x) asinh(x), 
#             inverse = function(x) sinh(x))
# }

# print('Starting plot')
# df$Type <- factor(df$group_type, levels = c("E", "Spiral", "S0"))
# g2 <- ggplot(df,aes(x=lg_M_B,y=N_GC))+
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
# Pres<-summary(as.mcmc.list(jags.neg, vars="PRes"),quantiles=0.5)$quantiles
# Dispersion = sum(Pres^2)/(N-4)# beta.0, beta.1 , beta.2, and k, 4 parameters


# S.NB1<-ggs(jagssamples.nb ,family=c("beta"))
# S.NB2<-ggs(jagssamples.nb,family=c("size"))

# S.NB<-rbind(S.NB1,S.NB2,deparse.level=2)
# S.NB$Parameter<-revalue(S.NB$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1]),"beta.2"=expression(beta[2]),"size"="k"))


# g1<-ggs_density(S.NB)+
#   scale_colour_economist(guide="none")+
#   theme_hc()+
#   scale_fill_economist()+
#   theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
#   theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
#         legend.position="none",plot.title = element_text(hjust=0.5),
#         axis.title.y=element_text(vjust=0.75),
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15),
#         strip.text.x=element_text(size=15),
#         axis.title.x=element_text(vjust=-0.25),
#         text = element_text(size=25))+xlab("Parameter  value")+ylab("Density") +
#   facet_grid(Parameter~.,labeller=label_parsed,scales = "free")

# CairoPDF("Figures/posterior_full_nb2.pdf",height=10,width=8)
# #facet_wrap_labeller(g1,labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),"k"))
# g1
# dev.off()

# g0<-ggs_traceplot(S.NB)+
#   scale_colour_economist(guide="none")+
#   theme_hc()+
#   scale_fill_economist()+
#   #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
#   theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
#         legend.position="none",plot.title = element_text(hjust=0.5),
#         axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
#         strip.text.x=element_text(size=25),
#         axis.title.x=element_text(vjust=-0.25),
#         text = element_text(size=25))+
#   ylab("Parameter value")+
#   xlab("Iteration")+
#   facet_grid(Parameter~.,labeller=label_parsed,scales = "free")

# CairoPDF("Figures/chain_full_nb2.pdf",height=10,width=8)
# g0 
# dev.off()

# jags.DIC <- jags.model(
#   data = jags.data, 
#   inits = inits1, 
#   textConnection(model.NB2),
#   n.chains = 3,
#   n.adapt=2000
# )

# update(jags.DIC , 5000)
# dicsamples.nb <- dic.samples(jags.DIC, params, n.iter = 10000,type="pD")

# dicsamples.nb





