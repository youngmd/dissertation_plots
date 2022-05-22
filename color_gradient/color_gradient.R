library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
require(plyr)
parseGccs <- function(datapath, galaxy){
  cnames1 <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr","raddist")

  cnames2 <- c("phot_id","bx","by","amb","bi","bi_err","vx","vy","amv","vi","vi_err","rx","ry","amr","ri","ri_err","vmag","bmv","vmr")

  gccdata <- tryCatch({
      read.table(datapath, col.names=cnames1)
    },
    warning = function(w){
      print()
    },
    error = function(err){
      return(read.table(datapath, col.names=cnames2))
    }
  )
  gccdata$bmr <- gccdata$bmv + gccdata$vmr
  gccdata$galname <- galaxy
  return(data.frame(bmr = gccdata$bmr, raddist = gccdata$raddist, galname= gccdata$galname))
}

parseGccs_simple <- function(datapath, galaxy){
  cnames1 <- c("bmr","raddist")
  gccdata <- read.table(datapath, col.names=cnames1)
  gccdata$galname <- galaxy
  return(gccdata)
}


lm_eqn <- function(df){
    m = lm(y ~ x, df);
    print(summary(m))
    eq <- substitute(italic(slope) == a~"±"~b, 
         list(a = format(coef(m)[2], digits = 2), 
              b = format(summary(m)$coefficients['x',"Std. Error"], digits = 2)))
    as.character(as.expression(eq));                 
}

plotColorGradient <- function(gccs){
  gccs <- subset(gccs, bmr < 1.9)
  gccs$x <- gccs$raddist
  gccs$y <- gccs$bmr
  eq <- ddply(gccs,.(galname),lm_eqn)
  g <- ggplot(gccs, aes(raddist,bmr)) + 
    geom_point(size=1, shape=22, alpha=0.2, color="black") +
    geom_smooth(method=lm) +
    geom_text(data=eq,aes(x = 1, y = 1.6,label=V1), parse = TRUE, inherit.aes=FALSE, hjust=0) +
    facet_wrap(~galname,nrow=7) +
    theme_bw() + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    # annotate("text", x = Inf, y = Inf, label = galname, size = 5, family="CM Roman",hjust = 1, vjust = 1) +
    ylab("B - R") + 
    xlab("Projected Radial Distance (arcmin)")
}

galname <- "NGC 5846"
datapath <- "~/gcc/n5846/gcfinder2/gc_cand_calib.90.out"
gccs <- parseGccs(datapath, galname)
#g1 <- plotColorGradient(gccs, galname)


galname <- "NGC 4649"
datapath <- "~/gcc/n4649/gcfinder2.n4649/gc_cand_calib.90.out"
gccs <- rbind(gccs, parseGccs(datapath, galname))
#g2 <- plotColorGradient(gccs, galname)

galname <- "NGC 4621"
datapath <- "~/gcc/n4649/gcfinder2.n4621/gc_cand_calib.90.out"
gccs <- rbind(gccs, parseGccs(datapath, galname))
#g3 <- plotColorGradient(gccs, galname)

galname <- "NGC 4382"
datapath <- "~/gcc/n4382/gcfinder2/gc_cand_calib.90.out"
gccs <- rbind(gccs, parseGccs(datapath, galname))
#g4 <- plotColorGradient(gccs, galname)

galname <- "NGC 1023"
datapath <- "~/gcc/n1023/gc_cand_calib.90.out"
gccs <- rbind(gccs, parseGccs(datapath, galname))

galname <- "NGC 7332"
datapath <- "~/gcc/n7332/gcc90.dat"
gccs <- rbind(gccs, parseGccs_simple(datapath, galname))

galname <- "NGC 4013"
datapath <- "~/gcc/n4013/gc_cand_calib.after_cut.90_percent_sample"
gccs <- rbind(gccs, parseGccs(datapath, galname))


#output to pdf
gplot <- plotColorGradient(gccs)
ggsave("color_gradient.pdf", gplot, width=8, height=8)
