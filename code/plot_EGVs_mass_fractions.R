#-----script header---------------------------------------------------------####
# Date: 25.10.2021
# Author: Diana Szeliova
# Description: plot mass fractions for each EGV
# Input: "concentrations_LDs_min_lip.csv" - metabolite concentrations [g/g]
# Output: "concentrations_EGVs_LDs_min_lip.png" 

#-----library declarations--------------------------------------------------####
library(rstudioapi)
library(RColorBrewer)

#-----read and process data-------------------------------------------------####

# set working directory to where the code is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder <- "../plots/"

conc <- read.csv("../data/concentrations_LDs_min_lip.csv", row.names = 1)


#----Parameters-------------------------------------------------------------####
xlim <- c(0, max(conc$mu))
lwd <- 1.6
oma <- c(2,5,1,0.5)
fig_size <- c(26,13)
leg_size <- 1.4

new_names_conc <- c()
for(n in 1:11){
  letter <- ifelse(n < 6, "v", "w")
  name_rxn <- gsub("R_", "", colnames(data)[n])
  new_names_conc <- c(new_names_conc, as.expression(bquote('x'[.(colnames(conc)[n])])))
}

colors <- data.frame(colors = c(brewer.pal(8, "Paired"), brewer.pal(3, "Dark2")),
                     mets = colnames(conc)[1:11],
                     lty = c(rep(1,5), rep(2,6)))


#----Functions--------------------------------------------------------------####
add_axes <- function(xlim, xtick, ylim, ytick, plot_name){
  x_axes <- c(paste("EGV", 25:26, sep = ""))
  y_axes <- c(paste("EGV", 0:4*8+1, sep = ""))
  size <- 1.2
  if(plot_name %in% x_axes){
    axis(1, seq(0, xlim[2], xtick), cex.axis = size)
  }
  if(plot_name %in% y_axes){
    axis(2, seq(0, ylim[2], ytick), cex.axis = size)
  }
  box()
}


add_labels_and_legend <- function(colors, ylab, legend_x, legend_y, lty, 
                                  new_names){
  title(xlab = bquote("Growth rate [h"^"-1"*"]" ), ylab = NA, outer = TRUE, 
        line = -6.5, cex.lab = 1.9)
  title(xlab = NA, ylab = ylab, outer = TRUE, 
        line = 2.5, cex.lab = 1.9)
  
  for(l in c(1:length(new_names))){
    legend(l*1.2-legend_x, legend_y, legend = new_names[(l+l-1):(l+l)],
           col = colors$colors[(l+l-1):(l+l)],
           lty = colors$lty[(l+l-1):(l+l)], 
           lwd = lwd, bty = "n", xpd = NA, cex = 1.8)
  }
}


#----Plot concentrations----------------------------------------------------####
png(filename = paste(folder, "concentrations_EGVs_LDs_min_lip.png", sep = ""), 
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)

ylim <- c(0, 1.15)
par(mfrow = c(4,8), oma = oma, mar = c(0,0,0,0))
for(egv in paste("EGV", 1:26, sep = "")){
  one_egv <- conc[conc$EGVs == egv,]
  for(column in colnames(conc)[1:11]){
    plot(one_egv$mu, one_egv[, column], lwd = lwd, type = "l",
         lty = colors[colors$mets == column, "lty"],
         xlim = xlim, ylim = ylim,
         xlab = NA, ylab = NA, axes = FALSE,
         col = colors[colors$mets == column, "colors"])
    par(new=TRUE)
  }
  box()
  letter <- ifelse(grepl("R_", column), "w", "v")
  add_axes(xlim, 0.5, ylim, 0.5, egv)
  legend(-0.3, ylim[2]*1.1, egv, 
         bty = "n", cex = leg_size, text.font = 2)
  abline(v=1.2131, col="grey46")
  par(new=FALSE)
}

add_labels_and_legend(colors, bquote("Mass fractions"),
                      -0.5, 0.64, 1, new_names_conc)
dev.off()