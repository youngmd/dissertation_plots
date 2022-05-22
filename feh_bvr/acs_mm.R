library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(mclust)
# Read in the data

options(warn=-1)

gcc <- read.table('acsvcs.dat',header=TRUE,colClasses=c("character", rep("numeric",3)))

gcc$gz<- gcc$gmag - gcc$zmag

BIC <- data.frame(galaxy='test', comp=0, BIC=0, type='E')
results <- data.frame(galaxy='test', n=0, G=0)
for(i in levels(factor(gcc$VCC))) { 

    gccdata <- subset(gcc, VCC == i)
    n <- nrow(gccdata)
    if(n < 400){
        next
    }
    #m <- Mclust(gccdata$gz)
    le <- mclustBootstrapLRT(gccdata$gz, modelName='E')
    bestG <- max(le$G)
    #lv <- mclustBootstrapLRT(gccdata$BR, modelName='V')
    df <- data.frame(galaxy=paste0("VCC",i), n=n, G=bestG)
    results <- rbind(results, df)
    }

results <- subset(results, galaxy != 'test')
