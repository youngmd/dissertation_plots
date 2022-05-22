# library(Cairo)
# library(ggplot2)
# library(ggthemes)
# library(scales)
# library(grid)
# library(extrafont)
# library(gtable)
# library(kSamples)
# 
# parseGccs <- function(datapath, galaxy){
#   cnames1 <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr","raddist")
# 
#   cnames2 <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr")
# 
#   gccdata <- tryCatch({
#       read.table(datapath, col.names=cnames1)
#     },
#     warning = function(w){
#       print()
#     },
#     error = function(err){
#       return(read.table(datapath, col.names=cnames2))
#     }
#   )
#   gccdata$bmr <- gccdata$bmv + gccdata$vmr
#   gccdata$galname <- galaxy
#   return(data.frame(bmr = gccdata$bmr, galname= gccdata$galname))
# }
# 
# parseGccs_simple <- function(datapath, galaxy){
#   cnames1 <- c("bmr","raddist")
#   gccdata <- read.table(datapath, col.names=cnames1)
#   gccdata$galname <- galaxy
#   return(data.frame(bmr = gccdata$bmr, galname= gccdata$galname))
# }
# 


gcc90<- read.csv('gcc90.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

gccdata <- subset(gcc90, galaxy_name == 'NGC5846')
galname <- "NGC 5846"
datapath <- "~/n5846.mm.pdf"
source('mm_plots2.R')

gccdata <- subset(gcc90, galaxy_name == 'NGC5813')
galname <- "NGC 5813"
datapath <- "~/n5813.mm.pdf"
source('mm_plots2.R')

gccdata <- subset(gcc90, galaxy_name == 'NGC4406')
galname <- "NGC 4406"
datapath <- "~/n4406.mm.pdf"
source('mm_plots2.R')

gccdata <- subset(gcc90, galaxy_name == 'NGC4594')
galname <- "NGC 4594"
datapath <- "~/n4594.mm.pdf"
source('mm_plots2.R')

gccdata <- subset(gcc90, galaxy_name == 'NGC4382')
galname <- "NGC 4382"
datapath <- "~/n4382.mm.pdf"
source('mm_plots2.R')

# #gccs <- parseGccs(datapath, galname)
# #g1 <- plotColorGradient(gccs, galname)
# 
# 
# galname <- "NGC 4649"
# datapath <- "~/gcc/n4649/gcfinder2.n4649/gc_cand_calib.90.out"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# #g2 <- plotColorGradient(gccs, galname)
# 
# galname <- "NGC 4621"
# datapath <- "~/gcc/n4649/gcfinder2.n4621/gc_cand_calib.90.out"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# #g3 <- plotColorGradient(gccs, galname)
# 
# galname <- "NGC 1172"
# datapath <- "~/gcc/Files_for_GCS_Database/n1172/color_sample90.inner_3.45"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 4406"
# datapath <- "~/gcc/Files_for_GCS_Database/n4406/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 4472"
# datapath <- "~/gcc/Files_for_GCS_Database/n4472/gc_cand.color_sample90"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 5813"
# datapath <- "~/gcc/Files_for_GCS_Database/n5813/NGC5813_Gcc90.tsv"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# galname <- "NGC 1023"
# datapath <- "~/gcc/n1023/gc_cand_calib.90.out"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 7332"
# datapath <- "~/gcc/n7332/gc_cand_calib.90.dat"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 4013"
# datapath <- "~/gcc/n4013/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 891"
# datapath <- "~/gcc/Files_for_GCS_Database/n891/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# galname <- "NGC 2683"
# datapath <- "~/gcc/Files_for_GCS_Database/n2683/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# galname <- "NGC 3379"
# datapath <- "~/gcc/Files_for_GCS_Database/n3379/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# galname <- "NGC 3384"
# datapath <- "~/gcc/Files_for_GCS_Database/n3384/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# galname <- "NGC 3556"
# datapath <- "~/gcc/Files_for_GCS_Database/n3556/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# galname <- "NGC 4157"
# datapath <- "~/gcc/Files_for_GCS_Database/n4157/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 4473"
# datapath <- "~/gcc/Files_for_GCS_Database/n4473/color_sample90.inner6.7"
# gccs <- rbind(gccs, parseGccs_simple(datapath, galname))
# 
# 
# galname <- "NGC 7331"
# datapath <- "~/gcc/Files_for_GCS_Database/n7331/gc_cand.color_sample90"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# galname <- "NGC 7814"
# datapath <- "~/gcc/Files_for_GCS_Database/n7814/gc_cand_calib.after_cut.90_percent_sample"
# gccs <- rbind(gccs, parseGccs(datapath, galname))
# 
# 
# gccs <- subset(gccs, bmr < 1.8)
# gccs <- subset(gccs, bmr > 0.9)
# 
# step = 0.02
# start = min(gccs$bmr)
# end = max(gccs$bmr)
# offset = (end - start) * 0.10
# start = start + offset
# end = end - offset
# 
# adtests <- data.frame(side=character(), cut=double(), adtest=double())
# 
# while (start < end)
# {
# 
#     lgccs <- subset(gccs, bmr < start)
#     rgccs <- subset(gccs, bmr > start)
# 
#     lgccs$galname <- factor(lgccs$galname)
#     rgccs$galname <- factor(rgccs$galname)
# 
#     splitl <- split(lgccs$bmr, f=lgccs$galname)
#     splitr <- split(rgccs$bmr, f=rgccs$galname)
# 
#     test <- ad.test(splitl)
#     newrow = data.frame(side='left', cut=start, adtest=test$ad[6])
#     adtests = rbind(adtests, newrow)
# 
#     test <- ad.test(splitr)
#     newrow = data.frame(side='right', cut=start, adtest=test$ad[6])
#     adtests = rbind(adtests, newrow)
#     start = start + step
# }
# 
# ad.test(split(gccs$bmr, f=gccs$galname))
# 
# ggplot(data=adtests, aes(x=cut, y=adtest, color=side, linetype=side)) + geom_line()
# 
# 
