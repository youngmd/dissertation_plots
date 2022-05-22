library(Cairo)
library(ggplot2)
library(extrafont)
# Read in the data
data = read.table('escut_results.dat', col.names=c("Filter","mag","fwhm","good"))

# Time out the bad stuff
gooddata <- subset(data, fwhm < 13)

# Set the order for facets
gooddata$Filter_s = factor(gooddata$Filter, levels=c('B','V','R'))

# create the plot
g <- ggplot(gooddata, aes(mag,fwhm)) + 
    geom_point(size=1, shape=22, alpha=0.2, color="grey40") +
    geom_point(data = subset(gooddata, good > 0), size=0.5, shape=21, alpha=0.2, color="black") +
    facet_wrap(~Filter_s,nrow=3) + 
    theme_bw() + 
    theme(text=element_text(family="CM Roman")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    ylab("Gaussian FWHM (pixels)") + 
    xlab("Instrumental Magnitude") + 
    xlim(-8, 0) + 
    ylim(2, 12)


#output to pdf
ggsave("escut.pdf", g, width=7, height=8)
embed_fonts("escut.pdf", outfile="escut_embed.pdf")
dev.off()
