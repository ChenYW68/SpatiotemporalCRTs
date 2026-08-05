rm(list = ls())

source("./LoadPackages/RDependPackages.R")

ate_file <- "./result/case/all_ATE_orginal_scale.RData"
required_ate_objects <- c(
  "JSTVC_spIE_10.ATE",
  "JSTVC_spIE_30.ATE",
  "JSTVC_spIE_50.ATE"
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

Title.1 <- TeX("(A) Spatial IEs (range = 10km)")
ATE.1 <- colSums(JSTVC_spIE_10.ATE)

Title.2 <- TeX("(B) Spatial IEs (range = 30km)")
ATE.2 <- colSums(JSTVC_spIE_30.ATE)

Title.3 <- TeX("(C) Spatial IEs (range = 50km)")
ATE.3 <- colSums(JSTVC_spIE_50.ATE)

ranking.1 <- order(as.numeric(ATE.1))
ranking.2 <- order(as.numeric(ATE.2))
ranking.3 <- order(as.numeric(ATE.3))

group_names <- group_names[order(ATE.1)]
color.arm <- data.frame(arms = group_names, color = Colo)

causal_effect.1 <- ATE.1[order(ATE.1)]
causal_effect.2 <- ATE.2[order(ATE.1)]
causal_effect.3 <- ATE.3[order(ATE.1)]

xlim.1 <- c(-1.0, 0)
xlim.2 <- c(-1.0, 0)
xlim.3 <- c(-1.0, 0)

text_shift.1 <- 0.07
text_shift.2 <- 0.07
text_shift.3 <- 0.07

panel.width <- c(0.75, 0.73, 0.70, 0.55)

draw_rank_panel <- function(title_text,
                            effects,
                            ranking,
                            x_limits,
                            text_shift,
                            group_names,
                            color_arm,
                            show_ranking_label = FALSE) {
  bar_height <- 0.6

  plot(
    NA,
    xlim = x_limits,
    ylim = c(0.5, 6),
    xlab = "",
    ylab = "",
    main = "",
    axes = FALSE
  )
  axis(1, labels = FALSE, tck = 0)
  mtext(title_text, side = 3, line = 2.6, adj = 0.1, cex = 1.3)
  mtext(expression(ATE == DE + IE), side = 3, line = 0.2, adj = 0.95, cex = 1)

  if (show_ranking_label) {
    mtext("Ranking of arms: ", side = 1, line = 1.1, adj = 0.03, cex = 1.5)
    ranking_adj <- 0.95
  } else {
    ranking_adj <- 0.72
  }

  mtext(
    paste0(
      "(",
      ranking[1], ", ", ranking[2], ", ", ranking[3], ", ",
      ranking[4], ", ", ranking[5], ", ", ranking[6], ")"
    ),
    side = 1,
    line = 1.1,
    adj = ranking_adj,
    cex = 1.5
  )

  for (i in seq_along(group_names)) {
    y_center <- i
    ind <- which(color.arm$arms %in% group_names[i])
    rect(
      xleft = effects[i],
      ybottom = y_center - bar_height / 2,
      xright = 0,
      ytop = y_center + bar_height / 2,
      col = color_arm$color[ind],
      border = NA
    )
    text(
      x = effects[i] + text_shift,
      y = y_center,
      labels = sprintf("%.2f", effects[i]),
      col = "black",
      adj = 0.5
    )
  }
}

pdf(
  file = file.path(figure_dir, "FigS16_Ranks_JSTVC_spIEs.pdf"),
  width = 15,
  height = 5
)

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

draw_rank_panel(
  title_text = Title.1,
  effects = causal_effect.1,
  ranking = ranking.1,
  x_limits = xlim.1,
  text_shift = text_shift.1,
  group_names = group_names,
  color_arm = color.arm,
  show_ranking_label = TRUE
)

draw_rank_panel(
  title_text = Title.2,
  effects = causal_effect.2,
  ranking = ranking.2,
  x_limits = xlim.2,
  text_shift = text_shift.2,
  group_names = group_names,
  color_arm = color.arm
)

draw_rank_panel(
  title_text = Title.3,
  effects = causal_effect.3,
  ranking = ranking.3,
  x_limits = xlim.3,
  text_shift = text_shift.3,
  group_names = group_names,
  color_arm = color.arm
)

par(mar = c(5, 0, 4, 0), mgp = c(3, 0.8, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.5, 3.5))
usr <- par("usr")
label_x <- usr[1] + 0.1
bar_height <- 0.6
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
