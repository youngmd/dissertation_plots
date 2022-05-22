library(Cairo)
library(ggplot2)
library(ggthemes)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(gtable)
library(bayesmix)
# Read in the data

options(warn=-1)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  library(extrafont)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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
bmr <- gccdata$bmv + gccdata$vmr


bayesmod <- BMMmodel(gccdata$bmr, k=1, priors=list(kind = "independence", parameter = "priorsUncertain", hierarchical = NULL)) # k=2 for two components

jcontrol <- JAGScontrol(variables = c("mu", "tau", "eta", "S"), burn.in = 1000, n.iter = 10000, seed = 10)

z <- JAGSrun(gccdata$bmr, model = bayesmod, control = jcontrol, tmp = FALSE, cleanup = TRUE)

zSort <- Sort(z, by = "mu")

zSort

