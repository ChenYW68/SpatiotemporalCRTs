rm(list = ls())
source("./LoadPackages/RDependPackages.R")
pdf(file  = "./figure/FigS11_Gneiting_monitor_convergence.pdf",
    width = 18,
    height = 7)
par(mar = c(3, 4, 2, 0), mgp = c(2, 0.8, 0))
par(mfrow = c(1, 2),
    cex = 1.5,
    cex.axis = 1.5,
    cex.lab = 1.5,
    cex.main = 1.5,
    lwd = 3)
load(paste0("./result/single_random/sim_random_MCMC_FALSE_250.RData"))
plot(CV_Ranalysis[[1]]$detailed.Para.List$ELBO/1e5,
     cex = 1,
     type = "l",
     # pch = 20,
     col = "black",
     xlab = "Iteration",
     ylab = expression(Evidence~Lower~Bound~"(×"~10^{5}*")"),
     xlim = c(1, 230),
     ylim = c(2, 4),
     axes = FALSE,
     frame.plot = FALSE)

# Add only the bottom and left axes
axis(1)  # x-axis
axis(2)  # y-axis
mtext("(A)", side = 3, line = 0.51, adj = -0.15, cex = 2)
load(paste0("./result/single_random/sim_random_MCMC_FALSE_298.RData"))
plot(CV_Ranalysis$detailed.Para.List$ELBO/1e5,
     cex = 1,
     type = "l",
     # pch = 20,
     col = "black",
     xlab = "Iteration",
     ylab = expression(Evidence~Lower~Bound~"(×"~10^{5}*")"),
     xlim = c(1, 230),
     ylim = c(2, max(CV_Ranalysis$detailed.Para.List$ELBO/1e5)),
     axes = FALSE,
     frame.plot = FALSE)

# Add only the bottom and left axes
axis(1)  # x-axis
axis(2)  # y-axis
mtext("(B)", side = 3, line = 0.51, adj = -0.18, cex = 2)
dev.off()
