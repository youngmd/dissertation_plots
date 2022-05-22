library(ggplot2)
gccdata <- read.csv('n5846_gccdata.csv', header=TRUE)
gclfdat <- read.table('n5846_gclf.dat', header=TRUE)

z <- subset(gccdata, galcen_ratio > 0.99)
contamdist <- hist(z$V, breaks='scott')
# data$vi <- abs(data$HRV - 1091)
# data$mn <- data$vi^2 * data$Rad
# ggplot(data=data, aes(x=Rad, y=HRV))+geom_point()+geom_point(data=gccs, aes(x=Rad, y=HRV), color='blue')+scale_y_continuous(limits=c(-1000,3000))

dev_a0 <- 3.69275
dev_a1 <- -1.86610 

contam <- 0.98 # 

extent <- 13.3


devPAlt <- function(r, v){
  da0 <- dev_a0
  da1 <- dev_a1
  c <- contam
  lsigma = da0 + da1*(r^0.25)
  sigma = 10^lsigma
  probs = vector(,length(r))
  for(j in 1:length(sigma)){
    for(i in 1:length(contamdist$breaks)){
       if(contamdist$breaks[i+1] > v[j]){
         fraccontam = contamdist$density[i]
         lwr = contamdist$breaks[i]
         upr = contamdist$breaks[i+1]
         nx <- x[x < upr & x > lwr]
         compl <- approx(gclfdat$Vmag, gclfdat$compl, v[j])
         print(v[j])
         print(compl$y)
         fracgcc <- ( length(nx)/ length(x) )
         siggcc <- sigma[j] * fracgcc * compl$y
         sigcontam <- c * fraccontam
         prob <- siggcc / (siggcc + sigcontam)
         probs[j] <- prob
         break
       }
    }
  }
  return(probs)
}

gccdata$P_GC <- devPAlt(gccdata$galcen_radius, gccdata$V)
