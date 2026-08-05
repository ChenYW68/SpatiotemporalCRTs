rm(list = ls())
source("./LoadPackages/RDependPackages.R")

pdf(file  = "./figure/FigS26_monitor_convergence_realData.pdf",
    width = 22,
    height = 6)

par(mar = c(3, 4, 2, 1), mgp = c(2, 0.8, 0))
par(mfrow = c(1, 3),
    cex = 1.4,
    cex.axis = 1.4,
    cex.lab = 1.4,
    cex.main = 1.4,
    lwd = 3)
# Kenya
load(paste0("./result/case/Kenya_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
plot(CV_Ranalysis$detailed.Para.List$ELBO/1e5,
     cex = 1,
     type = "l",
     # pch = 20,
     col = "black",
     xlab = "Iteration",
     ylab = expression(Evidence~Lower~Bound~"(×"~10^{5}*")"),
     xlim = c(1, 250),
     ylim = c(0, 0.8),
     axes = FALSE,
     frame.plot = FALSE)

axis(1)  # x-axis
axis(2)  # y-axis
mtext("(A)", side = 3, line = 0.51, adj = -0.18, cex = 2)
#Tanzania
load(paste0("./result/case/Tanzania_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
plot(CV_Ranalysis$detailed.Para.List$ELBO/1e5,
     cex = 1,
     type = "l",
     # pch = 20,
     col = "black",
     xlab = "Iteration",
     ylab = expression(Evidence~Lower~Bound~"(×"~10^{5}*")"),
     xlim = c(1, 250),
     ylim = c(0, 0.2),
     axes = FALSE,
     frame.plot = FALSE)

axis(1)  # x-axis
axis(2)  # y-axis
mtext("(B)", side = 3, line = 0.51, adj = -0.18, cex = 2)
# Both
load(paste0("./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
plot(CV_Ranalysis$detailed.Para.List$ELBO/1e5,
     cex = 1,
     type = "l",
     # pch = 20,
     col = "black",
     xlab = "Iteration",
     ylab = expression(Evidence~Lower~Bound~"(×"~10^{5}*")"),
     xlim = c(1, 250),
     ylim = c(0, 2),# max(CV_Ranalysis$detailed.Para.List$ELBO/1e5)
     axes = FALSE,
     frame.plot = FALSE)

axis(1)  # x-axis
axis(2)  # y-axis
mtext("(C)", side = 3, line = 0.51, adj = -0.18, cex = 2)

dev.off()
