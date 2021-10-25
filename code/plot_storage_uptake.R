#-----script header---------------------------------------------------------####
# Date: 19.10.2021
# Author: Diana Szeliova
# Description: plot carbon storage fraction of biomass (L+LD+G) and N uptake 
#             (v*omega/mu)
# Input: "EGVs_LDs_min_lip.csv" - fluxes [mmol/g*h]
#        "concentrations_LDs_min_lip.csv" - metabolite concentrations [g/g]
# Output: "storage_fraction_LDs_min_lip_groups_y.png"


#-----library declarations--------------------------------------------------####
library(rstudioapi)
library(RColorBrewer)


#-----read and process data-------------------------------------------------####
# set working directory to where the code is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder <- "../plots/"

data_l <- read.csv("../data/EGVs_LDs_min_lip.csv", row.names = 1)
conc_l <- read.csv("../data/concentrations_LDs_min_lip.csv", row.names = 1)

groups <- list(group1 = c("EGV1", "EGV2"),
               group2 = c("EGV3"),
               group3 = c(paste("EGV", 4:8, sep = "")),
               group4 = c("EGV9", "EGV10", "EGV17", "EGV18"),
               group5 = c("EGV11", "EGV19"),
               group6 = c(paste("EGV", 12:16, sep = ""), paste("EGV", 20:26, sep = "")))
graphical <- list(group1 = c("1a: EGV1-2", "blue", 3),
                group2 = c("2a: EGV3", "red", 3),
                group3 = c("3a: EGV4-8", "#00FFFF", 1),
                group4 = c("1b: EGV9-10/17-18", "blue", 4),
                group5 = c("2b: EGV11/19", "red", 4),
                group6 = c("3b: EGV12-16/20-24 + EGV25-26", "#00FFFF", 1))    


#----Parameters-------------------------------------------------------------####
xlim <- c(-0.05, max(data_l$mu)*1.05)
fig_size <- c(25,14)
ltys <- c(1,1:3)
lwd <- 2.5
axis_size <- 1.3

# ylims for top and bottom (zoomed in) parts of the plots
plot_parts <- list("high" = c(0.15,1.25), "low" = c(-0.01, 0.15))
panels <- list("upper" = 1:3, "lower" = 4:6)

#----Make plots-------------------------------------------------------------####

png(filename = paste(folder, "storage_fraction_LDs_min_lip_groups_y.png", sep = ""), 
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)
par(mfrow = c(4,3), mar = c(0,0,0,0), oma = c(3.5,4.5,0.5,0.5))

for(panel in panels){
  for(part in names(plot_parts))
    for(group in names(groups[panel])){
      for(egv in groups[[group]]){
        one_egv_c <- conc_l[conc_l$EGVs == egv,]
        one_egv <- data_l[data_l$EGVs == egv,]
        ratio <- one_egv$IN*0.018/one_egv$mu
        
        # adjust margins
        ylim <- plot_parts[[part]]
        if(part == "low"){
          par(mar = c(1,0,0.2,0))
        }else{
          par(mar = c(0.3,0,0.1,0))
        }
        
        # plot uptakes
        plot(one_egv$mu, ratio, axes = FALSE, yaxs="i", xaxs = "i",
             ylim = ylim, xlim = xlim, type = "l", lwd = 2.5, 
             lty = as.numeric(graphical[[group]][3]),
             col = graphical[[group]][2])
        
        par(new=TRUE)
        
        # select line type based on where carbon is stored
        if(max(one_egv_c$LD)<10e-16 & max(one_egv_c$G)<10e-16){
          lty <- ltys[2]
        }else if(max(one_egv_c$LD)>10e-16 & max(one_egv_c$G)<10e-16){
          lty <- ltys[3]
        }else if(max(one_egv_c$LD)<10e-16 & max(one_egv_c$G)>10e-16){
          lty <- ltys[3]
        }else{
          # nothing else should be there
          lty <- 4
        }
        
        # plot storage content
        plot(one_egv$mu, ((one_egv_c$L+one_egv_c$LD+one_egv_c$G)/rowSums(one_egv_c[1:11])), 
             col = "black", lty = lty, axes = FALSE, yaxs="i", xaxs = "i",
             ylim = ylim, xlim = xlim, type = "l", lwd = 2.5)
        
        # line that separates two phases
        abline(v=1.2131, col="grey69", lwd = 1.5, lty = 2)
        
        # add legend
        if(part == "high"){
          legend(-0.15, 1.3, graphical[[group]][1], bty = "n", 
                 cex = 1.4, text.font = 2, xjust=0)
        }

        # add lines
        if(part == "high" & egv %in% c(groups[["group1"]][1], groups[["group2"]][1], groups[["group3"]][1])){
          abline(h = ylim[1], lwd = lwd, col = "grey69")
          abline(h = ylim[2], lwd = lwd)
        }
        if(part == "low" & egv %in% c(groups[["group1"]][1], groups[["group2"]][1], groups[["group3"]][1])){
          abline(h = ylim[2], lwd = lwd, col = "grey69")
          abline(h = ylim[1], lwd = lwd)
        }
        if(part == "low" & egv %in% c(groups[["group4"]][1], groups[["group5"]][1], groups[["group6"]][1])){
          abline(h = ylim[2], lwd = lwd, col = "grey69")
          abline(h = ylim[1], lwd = lwd)
        }
        if(part == "high" & egv %in% c(groups[["group4"]][1], groups[["group5"]][1], groups[["group6"]][1])){
          abline(h = ylim[1], lwd = lwd, col = "grey69")
          abline(h = ylim[2], lwd = lwd)
        }
        abline(v = xlim[2], lwd = lwd)
        
        # add custom axes
        if(egv %in% c("EGV1", "EGV9")){
          abline(v=xlim[1], lwd = lwd)
          if(part == "high"){
            axis(2, at = seq(0, 1.5, 0.5), cex.axis = axis_size, lwd=0, lwd.ticks = 1)
          }else{
            axis(2, at = seq(-0.05, 0.15, 0.05), cex.axis = axis_size, lwd=0, lwd.ticks = 1)
          }
        }
        
        if(egv %in% c("EGV9", "EGV11", "EGV12") & part == "low"){
          axis(1, at = seq(-0.2, 1.4, 0.2), cex.axis = axis_size, lwd=0, lwd.ticks = 1)
        }
        
        par(new=TRUE)
      }
      par(new=FALSE)
  }
}


# outer axis titles
mtext(bquote("Growth rate [h"^"-1"*"]" ), side = 1, outer =  TRUE, 
      line = 2.6, cex = 1.25)
mtext(as.expression(bquote("Carbon storage / Ammonium uptake")), side = 2, outer =  TRUE, 
      line = 2.5, cex = 1.25)

dev.off()



