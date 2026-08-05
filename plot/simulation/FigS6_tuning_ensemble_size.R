load("./result/Simulation_300/random_tuning/VB_summary_Ne_Sensitibity.RData")
library(farver)
library(scatterplot3d)
sum_beta_full <- sum_beta[sum_beta$Ne >= 50 & sum_beta$Ne <= 500, ]
colors   <- colorRampPalette(c("blue", "red"))(10)
z_colors <- colors[cut(sum_beta_full$SD, breaks = 10)]
pdf(file  = "./figure/FigS6_Ne_Sensitivity.pdf",
    width = 16,
    height = 7)

par(mar = c(3.5, 4, 1, 1), mgp = c(2.5, 0.8, 0))
par(mfrow = c(1, 1),
    cex = 1.2,
    cex.axis = 1.3,
    cex.lab = 1.5,
    cex.main = 1.5,
    lwd = 1)
plot(
  sum_beta_full$Ne, sum_beta_full$RMSE,
  pch = 19, col = scales::alpha("blue", alpha = 0.5),
  xlab = expression(Ensemble~size~"(Ne)"),
  ylab = "Average MSE across regression coefficients", #expression(MSE~"(×"~10^{-6}*")")#,
)
dev.off()
