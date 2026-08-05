rm(list = ls())

source("./LoadPackages/RDependPackages.R")

ate_file <- "./result/case/all_ATE_orginal_scale.RData"
required_ate_objects <- c(
  "JSTVC_IE.ATE",
  "Sub_AR1.ATE",
  "Arm_AR1.ATE"
)

if (file.exists(ate_file)) {
  load(ate_file)
}

if (!all(required_ate_objects %in% ls())) {
  source("./plot/case/Ranking/0_all_ATE_orginal_scale.R")
}

figure_dir <- "./figure/"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

colour <- c("#2ca02c", "#1f77b4", "skyblue", "red", "#8856a7", "#ff7f0e")
c.alpha <- 0.3
Colo <- c(
  scales::alpha(colour[1], alpha = c.alpha),
  scales::alpha(colour[2], alpha = c.alpha),
  scales::alpha(colour[3], alpha = c.alpha),
  scales::alpha(colour[4], alpha = c.alpha),
  scales::alpha(colour[5], alpha = c.alpha),
  scales::alpha(colour[6], alpha = c.alpha)
)

group_names <- c(
  "Arm 1: CWT, CWT, CWT, CWT",
  "Arm 2: CWT, CWT, SBT, SBT",
  "Arm 3: CWT, CWT, No, No",
  "Arm 4: SBT, SBT, SBT, SBT",
  "Arm 5: SBT, SBT, No, No",
  "Arm 6: SBT, No, SBT, No"
)
color.arm <- data.frame(arms = group_names, color = Colo)

Title.1 <- TeX(" (A) JSTVC$_{-IE}$ ")
ATE.1 <- colSums(JSTVC_IE.ATE)

Title.2 <- TeX("(B) Sub-AR$_1$")
ATE.2 <- colSums(Sub_AR1.ATE)

Title.3 <- TeX("(C) Arm-AR$_1$")
ATE.3 <- colSums(Arm_AR1.ATE)

ranking.1 <- order(as.numeric(ATE.1))
ranking.2 <- order(as.numeric(ATE.2))
ranking.3 <- order(as.numeric(ATE.3))

group_names <- group_names[order(ATE.1)]
causal_effect.1 <- ATE.1[order(ATE.1)]
causal_effect.2 <- ATE.2[order(ATE.1)]
causal_effect.3 <- ATE.3[order(ATE.1)]

xlim.1 <- c(-0.3, 0.2)
xlim.2 <- c(-0.7, 0)
xlim.3 <- c(-0.7, 0)

text_shift.1.neg <- 0.025
text_shift.1.pos <- -0.029
text_shift.2 <- 0.045
text_shift.3.neg <- 0.045
text_shift.3.pos <- -0.06

panel.width <- c(0.55, 0.55, 0.55, 0.4)

pdf(file = file.path(figure_dir, "FigS13_Ranks_competing_methods.pdf"), width = 15, height = 5)

layout(
  matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4),
  widths = panel.width
)

par(
  mar = c(2.8, 0, 4, 0),
  mgp = c(3, 0.8, 0),
  cex = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.1,
  cex.main = 1,
  lwd = 1
)

plot(
  NA,
  xlim = xlim.1,
  ylim = c(0.5, 6),
  xlab = "",
  ylab = "",
  main = "",
  axes = FALSE
)
axis(1, labels = FALSE, tck = 0)
mtext(Title.1, side = 3, line = 2.6, adj = 0.1, cex = 1.3)
mtext(expression(ATE == DE), side = 3, line = 0.2, adj = 0.77, cex = 1)
mtext("Ranking of arms: ", side = 1, line = 1.1, adj = 0.05, cex = 1.5)
mtext(
  paste0(
    "(", ranking.1[1], ", ", ranking.1[2], ", ", ranking.1[3], ", ",
    ranking.1[4], ", ", ranking.1[5], ", ", ranking.1[6], ")"
  ),
  side = 1,
  line = 1.1,
  adj = 0.9,
  cex = 1.5
)

bar_height <- 0.6
for (i in seq_along(group_names)) {
  y_center <- i
  ind <- which(color.arm$arms %in% group_names[i])
  rect(
    xleft = causal_effect.1[i],
    ybottom = y_center - bar_height / 2,
    xright = 0,
    ytop = y_center + bar_height / 2,
    col = color.arm$color[ind],
    border = NA
  )
  text(
    x = ifelse(causal_effect.1[i] < 0, causal_effect.1[i] + text_shift.1.neg, causal_effect.1[i] + text_shift.1.pos),
    y = y_center,
    labels = sprintf("%.2f", causal_effect.1[i]),
    col = "black",
    adj = 0.4
  )
}

plot(
  NA,
  xlim = xlim.2,
  ylim = c(0.5, 6),
  xlab = "",
  ylab = "",
  main = "",
  axes = FALSE
)
axis(1, labels = FALSE, tck = 0)
mtext(Title.2, side = 3, line = 2.6, adj = 0.1, cex = 1.3)
mtext(expression(ATE == DE + IE), side = 3, line = 0.2, adj = 0.95, cex = 1)
mtext(
  paste0(
    "(", ranking.2[1], ", ", ranking.2[2], ", ", ranking.2[3], ", ",
    ranking.2[4], ", ", ranking.2[5], ", ", ranking.2[6], ")"
  ),
  side = 1,
  line = 1.1,
  adj = 0.5,
  cex = 1.5
)

for (i in seq_along(group_names)) {
  y_center <- i
  ind <- which(color.arm$arms %in% group_names[i])
  rect(
    xleft = causal_effect.2[i],
    ybottom = y_center - bar_height / 2,
    xright = 0,
    ytop = y_center + bar_height / 2,
    col = color.arm$color[ind],
    border = NA
  )
  text(
    x = causal_effect.2[i] + text_shift.2,
    y = y_center,
    labels = sprintf("%.2f", causal_effect.2[i]),
    col = "black",
    adj = 0.5
  )
}

plot(
  NA,
  xlim = xlim.3,
  ylim = c(0.5, 6),
  xlab = "",
  ylab = "",
  main = "",
  axes = FALSE
)
axis(1, labels = FALSE, tck = 0)
mtext(Title.3, side = 3, line = 2.6, adj = 0.1, cex = 1.3)
mtext(expression(ATE == DE + IE), side = 3, line = 0.2, adj = 0.95, cex = 1)
mtext(
  paste0(
    "(", ranking.3[1], ", ", ranking.3[2], ", ", ranking.3[3], ", ",
    ranking.3[4], ", ", ranking.3[5], ", ", ranking.3[6], ")"
  ),
  side = 1,
  line = 1.1,
  adj = 0.65,
  cex = 1.5
)

for (i in seq_along(group_names)) {
  y_center <- i
  ind <- which(color.arm$arms %in% group_names[i])
  rect(
    xleft = causal_effect.3[i],
    ybottom = y_center - bar_height / 2,
    xright = 0,
    ytop = y_center + bar_height / 2,
    col = color.arm$color[ind],
    border = NA
  )
  text(
    x = ifelse(causal_effect.3[i] < 0, causal_effect.3[i] + text_shift.3.neg, causal_effect.3[i] + text_shift.3.pos),
    y = y_center,
    labels = sprintf("%.2f", causal_effect.3[i]),
    col = "black",
    adj = 0.5
  )
}

par(mar = c(5, 0, 4, 0), mgp = c(3, 0.8, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.5, 3.5))
usr <- par("usr")
label_x <- usr[1] + 0.1
y.x <- c(-0.05, 0.58, 1.25, 1.88, 2.52, 3.2)

for (i in seq_len(6)) {
  ind <- which(color.arm$arms %in% group_names[i])
  text(
    x = label_x,
    y = bar_height / 2 + y.x[i],
    labels = group_names[i],
    adj = 0,
    col = scales::alpha(color.arm$color[ind], alpha = 0.5),
    xpd = TRUE,
    cex = 1
  )
}

dev.off()
