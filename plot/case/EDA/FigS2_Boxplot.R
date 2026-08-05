rm(list = ls())
#-----------------------------------------
source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site
Ken.G    <- G.mat

load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Tan.G    <- G.mat
Site     <- rbind(Ken.Site, Tan.Site)
#-----------------------------------------
region_flags <- c("Northwestern", "Northeastern", "Southern", "Western", "Eastern")
Ken_indices  <- lapply(region_flags[1:3], function(f) which(Ken.Site$flag == f))
Tan_indices  <- lapply(region_flags[4:5], function(f) which(Tan.Site$flag == f))
Score_Data   <- temp <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
Score_Data$Prevalence <- Score_Data$Prevalence*100

regions <- c("Northwestern", "Northeastern", "Southern", "Eastern", "Western")

# Dentification of outliers
df   <- subset(Score_Data, flag == regions[2])
temp <- df[df$Study_Arm %in% 1, c(1, 5, 7, 6, 8, 3)]
setorder(temp, -"Prevalence")
temp[1:5,]
# KEN212, KEN199



df   <- subset(Score_Data, flag == regions[3])
temp <- df[df$Study_Arm %in% 4, c(1, 5, 7, 6, 8, 3)]
setorder(temp, -"Prevalence")
temp[1:5,]
# KEN086


lab <- c("(A)", "(B)", "(C)", "(D)", "(E)")

pdf(file  = "./figure/FigS2_boxplots.pdf",
    width = 14,
    height = 10)


par(mfrow = c(2, 3),
    mar = c(3.5, 3.5, 2, 0.5),
    mgp = c(2, 0.8, 0),
    cex = 1.2,
    cex.axis = 1.3,
    cex.lab = 1.2,
    cex.main = 1,
    lwd = 1)






for (i in seq_along(regions)) {
  sa <- regions[i]
  df <- subset(Score_Data, flag == sa)
  bp <- boxplot(Prevalence ~ Study_Arm,
                data = df,
                ylab = "Prevalence (%)",
                xlab = "Arms",
                col = "lightblue",
                outline = T,
                ylim = c(0, 100)
  )
  mtext(lab[i], side = 3, line = 0.51, adj = 0, cex = 2)
}

dev.off()
