#-----script header---------------------------------------------------------####
# Date: 25.10.2021
# Author: Diana Szeliova
# Description: plot fluxes and mass fractions, grouped by reactions/metabolites
# Input: "EGVs_LDs_min_lip.csv" - fluxes [mmol/g*h]
#        "concentrations_LDs_min_lip.csv" - metabolite concentrations [g/g]
# Output: "fluxes_rxns_LDs_min_lip_cons.png"
#         "fluxes_rxns_LDs_min_lip_rest.png"
#         "concentrations_mets_LDs_min_lip_const.png"
#         "concentrations_mets_LDs_min_lip_rest.png"


#-----library declarations--------------------------------------------------####
library(rstudioapi)
library(RColorBrewer)

#-----read and process data-------------------------------------------------####

# set working directory to where the code is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder <- "../plots/"

data <- read.csv("../data/EGVs_LDs_min_lip.csv", row.names = 1)
conc <- read.csv("../data/concentrations_LDs_min_lip.csv", row.names = 1)

# convert fluxes from mmol/gh to g/gh
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

# split EGVs -- 8 EGVs that exist at all mus, 16/2 EGVs that exist at low/high mus 
consistent_modes <- c("EGV1", "EGV2", "EGV3", "EGV4", "EGV5", "EGV6", "EGV7", "EGV8")
data8 <- data[data$EGVs %in% consistent_modes,]
conc8 <- conc[conc$EGVs %in% consistent_modes,]
data_rest <- data[!data$EGVs %in% consistent_modes,]
conc_rest <- conc[!conc$EGVs %in% consistent_modes,]

#----Parameters-------------------------------------------------------------####
xlim <- c(0, max(data$mu))
oma <- c(4,4,1.5,0.5)
fig_size <- c(26,13)
leg_size <- 1.6
main_size <- 1.7
mu_crit <- 1.2131
lwd <- 1.5

set.seed(42)
EGVs_color_order <- c("EGV21", "EGV17", "EGV22", "EGV4", "EGV3", "EGV12", 
                      "EGV24", "EGV2", "EGV23", "EGV16", "EGV18", "EGV20", 
                      "EGV13", "EGV11", "EGV14", "EGV19", "EGV7", "EGV6", 
                      "EGV15", "EGV9", "EGV10", "EGV5", "EGV1", "EGV8", 
                      "EGV25", "EGV26")
colors <- data.frame(colors = sample(c(brewer.pal(10, "Paired"), brewer.pal(10, "Spectral"), 
                                brewer.pal(8, "Dark2")), 26),
                     EGVs = EGVs_color_order,
                     pch = c(rep(21,12), rep(22,12), rep(23,2)),
                     lty = c(rep(1,10), rep(2,8), rep(3,8)))
colors <- colors[match(paste("EGV", 1:26, sep = ""), colors$EGVs),]

#----Functions--------------------------------------------------------------####
make_plot <- function(data, column, colors, EGV, xlim, ylim){
  plot(data[, "mu"], data[, column], 
       lwd = 1.6, type = "l",
       xlim = xlim, ylim = ylim,
       xlab = NA, ylab = NA,
       cex.axis = 1.2,
       lty = colors[colors$EGVs == EGV, "lty"],
       col = colors[colors$EGVs == EGV, "colors"]) 
}


add_axes <- function(data, plot_name, x_axes, y_axes){
  size <- 1.2
  if(plot_name %in% x_axes){
    axis(1, seq(0, xlim[2], 0.2), cex.axis = size)
  }
  if(plot_name %in% y_axes){
    axis(2, seq(0, 1, 0.2), cex.axis = size)
  }
  box()
}


add_labels <- function(ylab){
  title(xlab = bquote("Growth rate [h"^"-1"*"]" ), ylab = NA, outer = TRUE, 
        line = 2.5, cex.lab = 2)
  title(xlab = NA, ylab = ylab, outer = TRUE, 
        line = 1.5, cex.lab = 2)
}


#----Plot fluxes of EGVs that exist at all mus------------------------------####
png(filename = paste(folder, "fluxes_rxns_LDs_min_lip_cons.png", sep = ""), 
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)

par(mfrow = c(3,4), oma = oma, mar = c(1.2,1.2,1.2,1.2))
for(column in colnames(data8)[1:11]){
  ylim <- c(0, max(data[, column]))
  for(EGV in unique(data8$EGVs)){
    one_EGV <- data8[data8$EGVs == EGV,]
    make_plot(one_EGV, column, colors, EGV, xlim, ylim)
    par(new=TRUE)
  }
  box()
  letter <- ifelse(grepl("R_", column), "w", "v")
  legend(-0.2, ylim[2]*1.1, as.expression(bquote(.(letter)[.(gsub("R_", "", column))])), 
         bty = "n", cex = leg_size, text.font = 2)
  abline(v = mu_crit, col = "grey55")
  par(new=FALSE)
}

add_labels(bquote("Fluxes [g g"^"-1" ~ "h"^"-1" * "]" ))

legend(1.45, 1.2, 
       legend = colors$EGVs[1:8],
       col = colors$colors[1:8],
       lty = colors$lty[1:8],
       lwd = lwd, bty = "n", xpd = NA, cex = 1.3)

title(main = "EGV1-8", cex.main = main_size, outer = TRUE)
dev.off()



#----Plot fluxes of EGVs that do not exist at all mus-----------------------####
png(filename = paste(folder, "fluxes_rxns_LDs_min_lip_rest.png", sep = ""), 
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)

par(mfrow = c(3,4), oma = oma, mar = c(1.2,1.2,1.2,1.2))
for(column in colnames(data_rest)[1:11]){
  ylim <- c(0,max(data[, column]))
  for(EGV in unique(data_rest$EGVs)){
    one_EGV <- data_rest[data_rest$EGVs == EGV,]
    make_plot(one_EGV, column, colors, EGV, xlim, ylim)
    par(new=TRUE)
  }
  box()
  letter <- ifelse(grepl("R_", column), "w", "v")
  legend(-0.2, ylim[2]*1.1, as.expression(bquote(.(letter)[.(gsub("R_", "", column))])), 
         bty = "n", cex = leg_size, text.font = 2)
  abline(v = mu_crit, col = "grey55")
  par(new=FALSE)
}

add_labels(bquote("Fluxes [g g"^"-1" ~ "h"^"-1" * "]" ))

for(l in c(1:2)){
  legend(l/1.5+0.78, 1.2, 
         legend = colors$EGVs[9:26][(l+8*(l-1)):((l+8*(l-1))+8)],
         col = colors$colors[9:26][(l+8*(l-1)):((l+8*(l-1))+8)],
         lty = colors$lty[9:26][(l+8*(l-1)):((l+8*(l-1))+8)],
         lwd = lwd, bty = "n", xpd = NA, cex = 1.3)
}

title(main = "EGV9-24 + EGV25-26", cex.main = main_size, outer = TRUE)
dev.off()



#----Plot concentrations that exist at all growth rates--------------------####
png(filename = paste(folder, "concentrations_mets_LDs_min_lip_const.png", sep = ""),
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)

par(mfrow = c(3,4), oma = oma, mar = c(1.2,1.2,1.2,1.2))
for(column in colnames(conc8)[1:11]){
  ylim <- c(0,max(conc[, column]))
  for(EGV in unique(conc8$EGVs)){
    one_EGV <- conc8[conc8$EGVs == EGV,]
    make_plot(one_EGV, column, colors, EGV, xlim, ylim)
    par(new=TRUE)
  }
  box()
  letter <- "x"
  legend("topleft", as.expression(bquote(.(letter)[.(column)])),
         bty = "n", cex = leg_size, text.font = 2)
  abline(v = mu_crit, col = "grey55")
  par(new=FALSE)
}
add_labels(bquote("Mass fractions"))

legend(1.47, 0.92, 
       legend = colors$EGVs[1:8],
       col = colors$colors[1:8],
       lty = colors$lty[1:8],
       lwd = lwd, bty = "n", xpd = NA, cex = 1.3)
title(main = "EGV1-8", cex.main = main_size, outer = TRUE)
dev.off()



#----Plot concentrations that do not exist at all growth rates---------------####
png(filename = paste(folder, "concentrations_mets_LDs_min_lip_rest.png", sep = ""),
    type="cairo", units="cm", width=fig_size[1], height=fig_size[2], res=300)

par(mfrow = c(3,4), oma = oma, mar = c(1.2,1.2,1.2,1.2))
for(column in colnames(conc_rest)[1:11]){
  ylim <- c(0,max(conc[, column]))
  for(EGV in unique(conc_rest$EGVs)){
    one_EGV <- conc_rest[conc_rest$EGVs == EGV,]
    make_plot(one_EGV, column, colors, EGV, xlim, ylim)
    par(new=TRUE)
  }
  box()
  letter <- "x"
  legend("topleft", as.expression(bquote(.(letter)[.(column)])),
         bty = "n", cex = leg_size, text.font = 2)
  abline(v = mu_crit, col = "grey55")
  par(new=FALSE)
}
add_labels(bquote("Mass fractions" ))

for(l in c(1:2)){
  legend(l/1.5+0.8, 0.92, 
         legend = colors$EGVs[9:26][(l+8*(l-1)):((l+8*(l-1))+8)],
         col = colors$colors[9:26][(l+8*(l-1)):((l+8*(l-1))+8)],
         lty = colors$lty[9:26][(l+8*(l-1)):((l+8*(l-1))+8)],
         lwd = lwd, bty = "n", xpd = NA, cex = 1.3)
}

title(main = "EGV9-24 + EGV25-26", cex.main = main_size, outer = TRUE)
dev.off()
