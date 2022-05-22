library(xtable)

gp <- read.csv('mastertable_gcs.csv', header=TRUE)

gp2 <- data.frame(Galaxy=gp$Galaxy)
rownames(gp2) <- NULL
gp2$NGC <- paste(gp$N_GC,"$\\pm$", gp$N_GC_err, sep=" ") 
gp2$SN <- paste(sprintf('%.1f', gp$S_N),"$\\pm$", sprintf('%.1f', gp$S_N_err), sep=" ") 
gp2$T <- paste(sprintf('%.1f', gp$T),"$\\pm$", sprintf('%.1f', gp$T_err), sep=" ") 
gp2$F_blue <- sprintf('%.2f', gp$F_blue)
gp2$Extent <- paste(sprintf('%.0f', gp$Extent),"$\\pm$", sprintf('%.0f', gp$Extent_err), sep=" ") 
gp2$Ref <- gp$Source

tmpTable <- xtable(gp2, caption ="Globular Clusters Stats", 
       label="table:mastergalprops") 

print(tmpTable, caption.placement="bottom", 
      sanitize.text.function= function(x){x}, include.rownames = FALSE) 