 rm(list = ls())
source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site
Ken.G    <- G.mat
load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Tan.G    <- G.mat
Site     <- rbind(Ken.Site, Tan.Site)
da  <- rbind(Kenya_Score_Data, Tanzania_Score_Data) %>%
    filter(Study_Arm %in% c(1, 2, 4),
           Year %in% c(2011:2015),
           !is.na(Prevalence))

pdf(file  = paste0("./figure/FigS4_Variance_intensity.pdf"),
    width = 12,
    height = 5)
par(mar = c(3.5, 3.5, 2, 0.5) + 0, mgp = c(2, 0.8, 0))
par(mfrow = c(1, 2),
    cex = 1.2,
    cex.axis = 1.3,
    cex.lab = 1.2,
    cex.main = 1,
    lwd = 1)

  mean.Intensity <- unique(da[, c(1, 40)])
  mean.Intensity[, 2] <- ((exp(mean.Intensity[, 2])))

  var.Intensity <- unique(da[, c(1, 41)])
  var.Intensity[, 2] <- ((exp(var.Intensity[, 2])))
  var.Prevalence <- aggregate(log(-log(1 - Prevalence)) ~ Village_ID,
                              data = da,
                              var, na.rm = T)

  colnames(var.Prevalence) <- c("Village_ID", "variance")
  var.Prevalence <- var.Prevalence %>%
    left_join(mean.Intensity, by = "Village_ID") %>%
    left_join(var.Intensity, by = "Village_ID")

  # var.Prevalence <- var.Prevalence[var.Prevalence$variance<3,]

  plot((var.Prevalence[, 3]), var.Prevalence[, 2], pch = 19,
       xlab = "Mean of log(intensity)", ylab = "Variance of log(-log(1 - prevalence))")
  mtext("(A)", side = 3, line = 0.51, adj = -0.17, cex = 2)
  plot((var.Prevalence[, 4]), var.Prevalence[, 2], pch = 19,
       xlab = "Variance of log(intensity)", ylab = "Variance of log(-log(1 - prevalence))")
  mtext("(B)", side = 3, line = 0.51, adj = -0.17, cex = 2)
dev.off()
























