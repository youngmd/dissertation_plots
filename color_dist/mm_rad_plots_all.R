library(ggplot2)
library(tidyr)
library(mclust)
library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(stargazer)
library(plyr)

# mwgc <-read.table('harris_clean.txt',header=TRUE)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc <- subset(mwgc, V.R != 0)
# mwgc$A_B <- mwgc$E.B.V. * 4.315
# mwgc$A_R <- mwgc$E.B.V. * 2.673
# mwgc$BR <- mwgc$B.V + mwgc$V.R - (mwgc$A_B - mwgc$A_R)
# mwgc <- subset(mwgc, BR < 1.6)

gcc90<- read.csv('gcc90.2.csv',header=TRUE)
gcc90<- subset(gcc90, BR < 1.8)
gcc90<- subset(gcc90, BR > 0.85)

gccall<- read.csv('gcc.all.csv',header=TRUE)
#gccall<- subset(gccall, X90percent != 1)

cnames <- c("r","N","sd","sderr","lsd","lsderrhi","lsderrlo","good_frac","r_d","NGC","Ncontam","contam_frac")
cnames2 <- c("r","N","sd","sderr","lsd","lsderrhi","lsderrlo","good_frac","r_d")


galpa = 30
galmajor = 6.9
galminor = 5.256103

nrp <- data.frame(class='A', galaxy='blah', x=0, sdr=0)

d <- subset(gcc90, galaxy_name == 'NGC4406')
d2 <- subset(gccall, galaxy_name == 'NGC4406')

ysize = 8619
d$bxt <- d$bx
d$bx <- d$by
d$by <- ysize - d$bxt
d2$bxt <- d2$bx
d2$bx <- d2$by
d2$by <- ysize - d2$bxt

extent = 17
galpa = 127.8
reff=1.56
galmajor = 11.482
galminor = 8.3176
datapath <- "~/gcc/Files_for_GCS_Database/n4406/"
rp = read.table(paste(datapath,'/radial_profile.dat',sep=''), col.names = cnames)
mask = read.table(paste0(datapath,'/rectregions.dat'), col.names=c('xcen','ycen','width','height'))
mask$x1 <- mask$xcen - 0.5*mask$width
mask$x2 <- mask$x1 + mask$width
mask$y1 <- mask$ycen - 0.5*mask$height
mask$y2 <- mask$ycen + mask$height
mask$factor <- as.character(row.names(mask))
mask$x1t <- mask$x1
mask$x1 <- mask$y1
mask$y1 <- ysize - mask$x1t
mask$x2t <- mask$x2
mask$x2 <- mask$y2
mask$y2 <- ysize - mask$x2t

outfile <- 'n4406_trimod.pdf'
galname <- 'NGC 4406 (E3)'
galcenx <- 4390
galceny <- 4370
bf <- 'dev'
d <- subset(d, galcen_radius < extent)
source('mm_rad_plots.R', print.eval=TRUE)
nrpt <- data.frame(class=newrp$class, galaxy=galname, x=newrp$r/extent, sdr=newrp$sdr)
nrp <- rbind(nrp,nrpt)
bf <- fits


d <- subset(gcc90, galaxy_name == 'NGC4594')
d2 <- subset(gccall, galaxy_name == 'NGC4594')

ysize = 8610
d$bxt <- d$bx
d$bx <- d$by
d$by <- ysize - d$bxt
d2$bxt <- d2$bx
d2$bx <- d2$by
d2$by <- ysize - d2$bxt


extent = 19
galpa = 89.5
reff=1.2
galmajor = 8.511
galminor = 4.8977
datapath <- "~/gcc/Files_for_GCS_Database/n4594/"
rp = read.table(paste(datapath,'/radial_profile.dat',sep=''), col.names = cnames)
mask = read.table(paste0(datapath,'/rectregions.dat'), col.names=c('xcen','ycen','width','height'))
mask$x1 <- mask$xcen - 0.5*mask$width
mask$x2 <- mask$x1 + mask$width
mask$y1 <- mask$ycen - 0.5*mask$height
mask$y2 <- mask$ycen + mask$height
mask$factor <- as.character(row.names(mask))

mask$x1t <- mask$x1
mask$x1 <- mask$y1
mask$y1 <- ysize - mask$x1t
mask$x2t <- mask$x2
mask$x2 <- mask$y2
mask$y2 <- ysize - mask$x2t

outfile <- 'n4594_trimod.pdf'
galname <- 'NGC 4594 (SAa)'
galcenx <- 4411
galceny <- 4375
d <- subset(d, galcen_radius < extent)
source('mm_rad_plots.R', print.eval=TRUE)
nrpt <- data.frame(class=newrp$class, galaxy=galname, x=newrp$r/extent, sdr=newrp$sdr)
nrp <- rbind(nrp,nrpt)
bf <- rbind(bf, fits)

d <- subset(gcc90, galaxy_name == 'NGC 5846')
d2 <- subset(gccall, galaxy_name == 'NGC 5846')
extent = 13.3
galpa = 47
reff=0.98
galmajor = 4.2657
galminor = 4.0742
datapath <- "~/gcc/n5846/gcfinder3"
rp = read.table(paste(datapath,'/corrected_profile.all.out',sep=''), col.names = cnames)
mask = read.table(paste0(datapath,'/rectregions.dat'), col.names=c('x1','x2','y1','y2'))
mask$factor <- as.character(row.names(mask))
outfile <- 'n5846_trimod.pdf'
galname <- 'NGC 5846 (E0-1)'
galcenx <- 3132
galceny <- 3896
d <- subset(d, galcen_radius < extent)
source('mm_rad_plots.R', print.eval=TRUE)
nrpt <- data.frame(class=newrp$class, galaxy=galname, x=newrp$r/extent, sdr=newrp$sdr)
nrp <- rbind(nrp,nrpt)
bf <- rbind(bf, fits)

d <- subset(gcc90, galaxy_name == 'NGC5813')
d2 <- subset(gccall, galaxy_name == 'NGC5813')
extent = 14.4
galpa = 143.3
reff=0.95
galmajor = 4.16
galminor = 2.6859
datapath <- "~/gcc/Files_for_GCS_Database/n5813"
rp = read.table(paste(datapath,'/corrected_profile.all.out',sep=''), col.names = cnames)
mask = read.table(paste0(datapath,'/rectregions.dat'), col.names=c('xcen','ycen','width','height'))
mask$x1 <- mask$xcen - 0.5*mask$width
mask$x2 <- mask$x1 + mask$width
mask$y1 <- mask$ycen - 0.5*mask$height
mask$y2 <- mask$ycen + mask$height
mask$factor <- as.character(row.names(mask))
outfile <- 'n5813_trimod.pdf'
galname <- 'NGC 5813 (E1-2)'
galcenx <- 3191
galceny <- 4214
d <- subset(d, galcen_radius < extent)
source('mm_rad_plots.R', print.eval=TRUE)
nrpt <- data.frame(class=newrp$class, galaxy=galname, x=newrp$r/extent, sdr=newrp$sdr)
nrp <- rbind(nrp,nrpt)
bf <- rbind(bf, fits)

d <- subset(gcc90, galaxy_name == 'NGC 4382')
d2 <- subset(gccall, galaxy_name == 'NGC 4382')
extent = 20.5
galpa = 30
reff=1.10
galmajor = 6.9
galminor = 5.256103
datapath <- "~/gcc/n4382/gcfinder3"
rp = read.table(paste(datapath,'/corrected_profile.all.out',sep=''), col.names = cnames)
mask = read.table(paste0(datapath,'/rectregions.dat'), col.names=c('x1','x2','y1','y2'))
mask$factor <- as.character(row.names(mask))
outfile <- 'n4382_trimod.pdf'
galname <- 'NGC 4382 (SA0pec)'
galcenx <- 3185
galceny <- 4620
d <- subset(d, galcen_radius < extent)
source('mm_rad_plots.R', print.eval=TRUE)
nrpt <- data.frame(class=newrp$class, galaxy=galname, x=newrp$r/extent, sdr=newrp$sdr)
nrp <- rbind(nrp,nrpt)
nrp <- subset(nrp, class != 'A')
bf <- rbind(bf, fits)

bf <- bf[c(8,1,2,3,5,4,6,7)]
stargazer(bf, summary=FALSE, rownames=FALSE,  digits = 2, column.labels = c("Galaxy", "Population", "Fit", "a0","a0err","a1","a1err","r_10"))