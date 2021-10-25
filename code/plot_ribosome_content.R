#-----script header---------------------------------------------------------####
# Date: 25.10.2021
# Author: Diana Szeliova
# Description: plot ribosome content in three groups
# Input:  "concentrations_LDs_min_lip.csv" - metabolite concentrations [g/g]
# Output: "ribosome_content.png"


#-----library declarations--------------------------------------------------####
library(rstudioapi)

#-----read and process data-------------------------------------------------####
# set working directory to where the code is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder <- "../plots/"

conc <- read.csv("../data/concentrations_LDs_min_lip.csv", row.names = 1)

groups <- list(group1 = c(paste("EGV", 1:4, sep = ""), paste("EGV", 9:12, sep = ""), paste("EGV", 17:20, sep = "")),
               group2 = c(paste("EGV", 5:7, sep = ""), paste("EGV", 13:15, sep = ""), paste("EGV", 21:23, sep = "")),
               group3 = c("EGV8", "EGV16", paste("EGV", 24:26, sep = "")))
               # group4 = c(paste("EGV", 9:12, sep = ""), paste("EGV", 17:20, sep = "")),
               # group5 = c(paste("EGV", 13:15, sep = ""), paste("EGV", 21:23, sep = "")),
               # group6 = c("EGV16", paste("EGV", 24:26, sep = "")))
graphical <- list(group1 = c("I: EGV1-4, 9-12, 17-20", "black", 1),
                  group2 = c("II: EGV5-7, 13-15, 21-23", "black", 1),
                  group3 = c("III: EGV8, 16, 24 + 25-26", "black", 1))
                  # group4 = c("II: EGV9-10/17-18", "blue", 4),
                  # group5 = c("III: EGV11/19", "red", 4),
                  # group6 = c("III: EGV12-16/20-24 + EGV25-26", "#00FFFF", 1))

#----Parameters-------------------------------------------------------------####
#colours <- c("black", brewer.pal(12, "Paired")[c(5)])
xlim <- c(-0.05, max(conc$mu))
ylim = c(0, 1)
fig_size <- c(25,14)
ltys <- c(1,1:3)
lwd <- 2.5
axis_size <- 1.3


png(filename = paste(folder, "ribosome_content.png", sep = ""), 
    type="cairo", units="cm", width=23, height=9, res=300)

par(mfrow = c(1,3), mar = c(0,0,0,0), oma = c(4.5,4,1,1))
for(group in names(groups)){
  for(egv in groups[[group]]){
    one_egv <- conc[conc$EGVs == egv,]
    
    plot(one_egv$mu, one_egv$R,
         xlim = xlim, ylim = ylim,
         type = "l", lty = 2, lwd = 1.5, col = "#7570B3",
         axes = FALSE)
    if(egv %in% c("EGV9", "EGV13", "EGV16")){
      axis(1, seq(0, 1.2, 0.2), cex.axis = 1.3)
    }
    if(egv %in% c("EGV1", "EGV9")){
      axis(2, seq(0, 1.2, 0.2), cex.axis = 1.3)
    }
    
    par(new=TRUE)
  }
  par(new=FALSE)
  
  abline(v=1.2131, col="grey69", lwd = 1.5, lty = 1)
  legend(-0.2, 1.05, graphical[[group]][1], bty = "n", 
         cex = 1.4, text.font = 2, xjust=0)
  box()
}


mtext(bquote("Growth rate [h"^"-1"*"]" ), side = 1, outer =  TRUE, 
      line = 3.3, cex = 1.25)
mtext(as.expression(bquote("Ribosome mass fraction")), side = 2, outer =  TRUE, 
      line = 2.4, cex = 1.25)

dev.off()
  
