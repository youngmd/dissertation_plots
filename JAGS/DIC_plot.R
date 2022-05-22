library(ggplot2)
library(ggthemes)
library(Cairo)

bod <- read.table('bod.dat', header=TRUE)
#bod$x = factor(bod$x2, levels=c('density','Mdyn','mass','R_e','L_K','L_V'))
bod$y = factor(bod$x1, levels=c('Mdyn','mass','R_e','L_K','L_V','density'))
bod$x = factor(bod$x2, levels=c('density','L_V','L_K','R_e','mass','Mdyn'))
bod$DIC <- floor(bod$DIC)

g <- ggplot(bod, aes(x,y)) + geom_tile(aes(fill=DIC), colour="white") + 
    scale_fill_gradient("DIC", high="white",low="steelblue2", limits=c(650,1200)) + 
    geom_text(aes(label= DIC)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)),
        axis.text.x=element_text(margin=margin(15,0,0,0)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.75, 0.75)) +
    scale_x_discrete(labels=c(expression(rho),expression(log~L[V]/L['\u0298']),expression(log~L[K]/L['\u0298']),expression(R[eff]),expression(log~M['*']/M['\u0298']),expression(log~M[dyn]/M['\u0298']))) + 
    scale_y_discrete(labels=c(expression(log~M[dyn]/M['\u0298']),expression(log~M['*']/M['\u0298']),expression(R[eff]),expression(log~L[K]/L['\u0298']),expression(log~L[V]/L['\u0298']),expression(rho)))
    #scale_x_discrete(labels=c(expression(rho),expression(log~M[dyn]/M['\u0298']),expression(log~M['*']/M['\u0298']),expression(R[eff]),expression(log~L[K]/L['\u0298']),expression(log~L[V]/L['\u0298'])))

