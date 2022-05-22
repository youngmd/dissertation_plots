library(ggplot2)

bod <- read.table('bod.dat', header=TRUE)

bod$x = factor(bod$x2, levels=c('L_K','L_V','Mdyn','mass','R_e'))
bod$y = factor(bod$x1, levels=c('R_e','mass','Mdyn','L_V','L_K'))

g <- ggplot(bod, aes(x, y))
g + geom_tile(aes(fill=DIC), colour="white") + 
    scale_fill_gradient(high="white",low="steelblue2") + 
    geom_text(aes(label=DIC)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=margin(0,15,0,0)),
        axis.text.x=element_text(margin=margin(15,0,0,0)),
        legend.position = c(0.75, 0.75)) +
    scale_y_discrete(labels=c(expression(R[eff]),expression(log~M['\u2726']/M['\u0298']),expression(log~M[dyn]/M['\u0298']),expression(log~L[V]/L['\u0298']),expression(log~L[K]/L['\u0298']))) +
    scale_x_discrete(labels=c(expression(log~L[K]/L['\u0298']),expression(log~L[V]/L['\u0298']),expression(log~M[dyn]/M['\u0298']),expression(log~M['\u2726']/M['\u0298']),expression(R[eff])))