#-----script header---------------------------------------------------------####
# Date:
# Author: Diana Szeliova
# Description: plot fluxes for each EGV
# Input: "EGVs_LDs_24_10.csv" - fluxes [mmol/g*h]
# Output: "fluxes_EGVs_LDs_min_lip_mw.png"

#-----library declarations--------------------------------------------------####
library(rstudioapi)
library(RColorBrewer)

#-----read and process data-------------------------------------------------####

# set working directory to where the code is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder <- "../plots/"

data <- read.csv("../data/EGVs_LDs_min_lip.csv", row.names = 1)


#----Parameters-------------------------------------------------------------####
xlim <- c(0, max(data$mu))
lwd <- 1.6
oma <- c(2,5,1,0.5)
fig_size <- c(26,13)
leg_size <- 1.4
cex <- 0.2

# masses
nL = 7   
nR = 7536*3  # multiplied by 3 to take into account the RNA mass of the ribosome
nI = 646
nE = 325
MW_Glc <- 0.18  # glucose MW [g/mmol]
MW_NH4 <- 0.018  # nitrogen MW [g/mmol]
MW_AA <- MW_Glc+MW_NH4
mw <- c(MW_Glc, MW_NH4, MW_AA, 
      nL*MW_Glc, nL*MW_Glc, 
      nI*MW_AA, nI*MW_AA, 
      nE*MW_AA, nE*MW_AA, nE*MW_AA, nR*MW_AA)

data[,1:11] <- data.frame(mapply(`*`,data[1:11], mw))

new_names_rxns <- c()
for(n in 1:11){
  letter <- ifelse(n < 6, "v", "w")
  name_rxn <- gsub("R_", "", colnames(data)[n])
  new_names_rxns <- c(new_names_rxns, as.expression(bquote(.(letter)[.(name_rxn)])))
}

colors <- data.frame(colors = c(brewer.pal(8, "Paired"), brewer.pal(3, "Dark2")),
                     rxns = colnames(data)[1:11],
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


add_labels_and_legend <- function(colors, ylab, legend_x, legend_y, 
                                  new_names){
  title(xlab = bquote("Growth rate [h"^"-1"*"]" ), ylab = NA, outer = TRUE, 
        line = -6.5, cex.lab = 1.9)
  title(xlab = NA, ylab = ylab, outer = TRUE, 
        line = 2.5, cex.lab = 1.9)
  
  for(l in c(1:length(new_names))){
    legend(l*1.3-legend_x, legend_y, legend = new_names[(l+l-1):(l+l)],
           col = colors$colors[(l+l-1):(l+l)],
           lty = colors$lty[(l+l-1):(l+l)], 
           lwd = lwd, bty = "n", xpd = NA, cex = 1.8)
  }
}


#----Plot fluxes------------------------------------------------------------####
png(filename = paste(folder, "fluxes_EGVs_LDs_min_lip_mw.png", sep = ""), 
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)

ylim <- c(0,max(data[,1:11]))
par(mfrow = c(4,8), oma = oma, mar = c(0,0,0,0))
for(egv in paste("EGV", 1:26, sep = "")){
  one_egv <- data[data$EGVs == egv,]
  for(column in colnames(data)[1:11]){
    plot(one_egv$mu, one_egv[, column], lwd = lwd, cex = cex, type = "l",
         lty = colors[colors$rxns == column, "lty"],
         xlim = xlim, ylim = ylim,
         xlab = NA, ylab = NA, axes = FALSE,
         col = colors[colors$rxns == column, "colors"])
    par(new=TRUE)
  }
  box()
  letter <- ifelse(grepl("R_", column), "w", "v")
  legend(-0.3, ylim[2]*1.1, egv, 
         bty = "n", cex = leg_size, text.font = 2)
  abline(v=1.2131, col="grey46")
  add_axes(xlim, 0.5, ylim, 0.5, egv)
  par(new=FALSE)
}


add_labels_and_legend(colors, bquote("Fluxes [g g"^"-1" ~ "h"^"-1" * "]" ),
                      -0.25, 0.7, new_names_rxns)
dev.off()


