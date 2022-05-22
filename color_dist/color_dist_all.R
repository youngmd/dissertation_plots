library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(gtable)
library(mixtools)
library(nortest)
library(ADGofTest)
library(fitdistrplus)
library(mclust)
# Read in the data

options(warn=-1)

lratioP <- function(crit, idf) {
  x <- crit
  an1 <- idf
  f <- x/an1
  an2 <- 1e10
  ff <- f
  prob <- 1.0

  if (f < 1.0) {
    ff <- 1.0 / f
    temp <- an1
    an1 <- an2
    an2 <- temp
  }

  a1 <- 2.0 / an1 / 9.0
  a2 <- 2.0 / an2 /9.0
  z <- abs(((1.0 - a2) * ff^0.3333 - 1.0 + a1)/sqrt(a2 * ff^0.666666 + a1))

  if (an2 <= 3.0) { z <- z*(1.0 + 0.08 * z^4 / an2^3)}

  fz <- exp(-z*z/2.0)*0.3989423
  w <-1.0 / (1.0 + z * 0.2316419)
  prob <- fz * w * ((((1.332074*w-1.821256)*w+1.781478)*w-0.3565638)*w + 0.3193815)

  if (f < 1.0) { prob <- 1.0 - prob}
  prob
}

gmm <- function(gccdata){

  mu2 = c(1.2, 1.4)
  l2 = c(0.5, 0.5)
  sig2 = c(0.1, 0.1)

  #m1 <- fitdist(gccdata$bmr, "norm")
  m2 <- normalmixEM(gccdata$BR, mu=mu2, maxit=10000, lambda=l2, sigma=sig2, k=2, mean.constr=c(NA, NA), arbmean=TRUE, arbvar=TRUE)

  ll1 <- logLik(glm(gccdata$BR~1))[1]

  m2$llr <- -2.0 * (ll1 - m2$loglik)
  m2$p <- lratioP(m2f_lr, 4)
  m2$n <- nrow(gccdata)

  return(m2)
}


gcc90<- read.csv('gcc90.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

# g1 <- ggplot(data=gccs, aes(BR)) + geom_histogram(binwidth = diff(range(gccs$BR))/20) + facet_wrap(~galaxy_name, scales = "free")

# results <- data.frame(galaxy='test', galtype='blah', galmass=0, mu1=0, mu2=0, l1=0, l2=0, sig1=0, sig2=0, llr=0, p=0, n=0)

# for(i in levels(gccs$galaxy_name)) { 
#     gccdata <- subset(gccs, galaxy_name == i)
#     if(nrow(gccdata) < 10){
#         print(paste("Skipping",i))
#         next
#     }
#     m <- gmm(gccdata)
#     l2 <- data.frame(galaxy=i, galtype=gccdata[1,]$host_galaxy_group_type, galmass=gccdata[1,]$host_galaxy_mass, mu1=m$mu[1], mu2=m$mu[2], l1=m$lambda[1], l2=m$lambda[2], sig1=m$sigma[1], sig2=m$sigma[2], llr=m$llr, p=m$p, n=m$n)
#     results <- rbind(results, l2)
# }

BIC <- data.frame(galaxy='test', comp=0, BIC=0, type='E')
results <- data.frame(galaxy='test', n=0, BIC='blah', LRT_E_n=0, LRT_E_p=0, LRT_V_n=0, LRT_V_p=0, ADGoF=0)
for(i in levels(gcc90$galaxy_name)) { 

    gccdata <- subset(gcc90, galaxy_name == i)
    n <- nrow(gccdata)
    if(n < 40){
        next
    }
    m <- Mclust(gccdata$BR)
    bestmod <- m$modelName
    bestG <- m$G
    validstr <- paste0(bestmod, bestG, ',')
    best <- m$bic
    ebics <- m$BIC[,'E']
    ebics <- subset(ebics, abs(ebics - best) < 5)
    vbics <- m$BIC[,'V']
    vbics <- subset(vbics, abs(vbics - best) < 5)
    for(j in names(ebics)){
        validstr <- paste0(validstr,"E",j)
    }
    for(j in names(vbics)){
        validstr <- paste0(validstr,"V",j)
    }
    #print(paste(i,validstr))
    # s <- data.frame(matrix(NA, nrow=nrow(m$BIC), ncol=4))
    # s <- rename(s, c("X1"='galaxy', 'X2'='comp', 'X3'='BIC', 'X4'='type'))
    # s$galaxy <- i
    # s$comp <- as.numeric(rownames(s))
    # s$BIC <- m$BIC[,1]
    # s$type <- "E"
    # BIC <- rbind(BIC, s)
    # s$BIC <- m$BIC[,2]
    # s$type <- "V"
    # BIC <- rbind(BIC, s)
    n <- nrow(gccdata)
    le <- mclustBootstrapLRT(gccdata$BR, modelName='E')
    lv <- mclustBootstrapLRT(gccdata$BR, modelName='V')
    df <- data.frame(galaxy=i, n=n, BIC=validstr, LRT_E_n=max(le$G), LRT_E_p=max(le$p.value), LRT_V_n=max(lv$G), LRT_V_p=max(lv$p.value), ADGoF=0)
    results <- rbind(results, df)
    }
results <- subset(results, galaxy != 'test')

# results <- data.frame(galaxy='test', BIC_n=0, BIC_m=0)
# for(i in levels(gcc90$galaxy_name)) { 
#     gccdata <- subset(gcc90, galaxy_name == i)
#     m <- Mclust(gccdata$BR)
#     if(m$G == 1){ next }
#     if(m$G > 5){ next}
#     df <- data.frame(galaxy=i, BIC_n=m$G, BIC_m=m$modelName)
#     results <- rbind(results, df)
#     }
# results <- subset(results, galaxy != 'test')
