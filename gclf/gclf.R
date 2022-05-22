library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(extrafont)
library(gtable)
# Read in the data

insert_minor <- function(major_labs, n_minor) {labs <- 
                              c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
                              labs[1:(length(labs)-n_minor)]}

cnames <- c("Vmag","oldnum","num","compl")

gclf = read.table(paste(datapath,'/corrected_gclf.out',sep=''), col.names=cnames)
bw <- gclf[2,]$Vmag-gclf[1,]$Vmag
gclf = subset(gclf, num>0)
count = sum(gclf$oldnum)

fitgclf = subset(gclf, compl > cutoff)
fitgclf = subset(fitgclf, oldnum > 3)

fitgclf$cnumerr <- sqrt(fitgclf$oldnum)/fitgclf$compl
xmin = mean - 5
xmax = mean + 5
dist <- data.frame( Vmag = seq(from = xmin, 
            to = xmax, 
            length.out = 500))

max_iters = 5000

all_chisq <- data.frame(norm = double(), chisq = double())


print("In here")
for (i in 1:max_iters) {
  fit <- dnorm(fitgclf$Vmag, mean=mean, sd=1.2) * i
  chisq = sum(((fitgclf$num - fit)^2)/(fitgclf$cnumerr)^2)
  newrow <- data.frame(norm=i,chisq=chisq)
  all_chisq = rbind(all_chisq, newrow)
}

a1 = all_chisq[which.min(all_chisq$chisq),]$norm

all_chisq <- data.frame(norm = double(), chisq = double())
for (i in 1:max_iters) {
  fit <- dnorm(fitgclf$Vmag, mean=mean, sd=1.3) * i
  chisq = sum(((fitgclf$num - fit)^2)/(fitgclf$cnumerr)^2)
  newrow <- data.frame(norm=i,chisq=chisq)
  all_chisq = rbind(all_chisq, newrow)
}


a2 = all_chisq[which.min(all_chisq$chisq),]$norm

all_chisq <- data.frame(norm = double(), chisq = double())
for (i in 1:max_iters) {
  fit <- dnorm(fitgclf$Vmag, mean=mean, sd=1.4) * i
  chisq = sum(((fitgclf$num - fit)^2)/(fitgclf$cnumerr)^2)
  newrow <- data.frame(norm=i,chisq=chisq)
  all_chisq = rbind(all_chisq, newrow)
}

a3 = all_chisq[which.min(all_chisq$chisq),]$norm

dist$dist1.2 <- dnorm(dist$Vmag,mean=mean, sd=1.2) * a1
dist$dist1.3 <- dnorm(dist$Vmag,mean=mean, sd=1.3) * a2
dist$dist1.4 <- dnorm(dist$Vmag,mean=mean, sd=1.4) * a3
dist$blank <- dnorm(dist$Vmag,mean=mean, sd=1.1) * a1 # for spacing at top of plot

if(max(dist$blank) < max(gclf$num)){
  dist$blank <- 1.1 * max(gclf$num)
}

peak = max(dist$dist1.2, 50)
offset = mean + 2
g <- ggplot() + 
    geom_bar(data=gclf, aes(x=Vmag, y=oldnum), stat="identity", width=bw)+
    geom_bar(data=gclf, aes(x=Vmag, y=num), stat="identity", width=bw, alpha=0.5)+
    # geom_histogram(data=gclf, aes(x=Vmag, weight=oldnum), binwidth=bw) +
    # geom_histogram(data=gclf, aes(x=Vmag, weight=num), binwidth=bw, alpha=0.5) +
    geom_path(data = dist, aes( Vmag, dist1.2), linetype=2, size=0.3, color='blue4') +
    geom_path(data = dist, aes( Vmag, dist1.3), linetype=3, size=0.3, color='green4') +
    geom_path(data = dist, aes( Vmag, dist1.4), linetype=4, size=0.3, color='red4') +
    geom_blank(data=dist, aes( Vmag, blank), stat="identity") + 
    theme_bw() + 
    theme(legend.position="none") + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-0.15, "cm"),
          axis.text.y=element_text(margin=margin(0,15,0,0)), 
          axis.text.x=element_text(margin=margin(15,0,0,0))) +
    ylab(expression('Number')) + 
    xlab(expression('V Magnitude')) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks=seq(0,40,1), labels=insert_minor(seq(0,40,2), 1)) + 
    annotate("text", x=Inf,y=Inf
, label = galname,size = 3,family="CM Roman", hjust = 1.2, vjust = 2)

#plot(mirror.ticks(g))

ggsave(outfile, g, width=7, height=7)
#embed_fonts('tmp.pdf', outfile=outfile)
