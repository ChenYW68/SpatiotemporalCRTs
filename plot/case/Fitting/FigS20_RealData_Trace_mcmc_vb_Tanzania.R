rm(list = ls())
source("./LoadPackages/RDependPackages.R")
load(paste0("./result/case/Tanzania_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
elapsed.vb.1 <- round(abs(run.time)[3]/3600, 1)
JSTVC.vb <- CV_Ranalysis$Process.monitoring$VB.EnKS.Iter#[, -2]
vb_long  <- JSTVC.vb %>%
  pivot_longer(
    cols = starts_with("pub."),
    names_to = "Coefficient",
    values_to = "Value"
  ) %>%
  mutate(Method = "VB")



load(paste0("./result/case/Tanzania_JSTVC_Xi_IEt_MCMC_pow_decay_loglog_11.RData")) #
elapsed.mcmc.1 <- round(abs(run.time)[3]/3600, 1)
JSTVC.mcmc <- CV_Ranalysis$Process.monitoring$VB.EnKS.Iter#[, -2]
mcmc_long  <- JSTVC.mcmc%>%
  filter(Iter >= 0e4) %>%
  pivot_longer(
    cols = starts_with("pub."),
    names_to = "Coefficient",
    values_to = "Value"
  ) %>%
  mutate(Method = "MCMC")

vb_long$type   <- "VB"
mcmc_long$type <- "MCMC"
trace_data <- bind_rows(vb_long, mcmc_long)

method_labels <- c(
  "MCMC" = paste0("MCMC~'(Tanzania; running time = ", elapsed.mcmc.1, "h)'"),
  "VB"   = paste0("VB~'(Tanzania; running time = ", elapsed.vb.1, "h)'")
)

coeff_labels <- c(
  "pub.Intercept" = "Intercept",
  "pub.CWT_1" = "CWT[1]",
  "pub.CWT_2" = "CWT[2]",
  "pub.CWT_3" = "CWT[3]",
  "pub.CWT_4" = "CWT[4]",
  "pub.SBT_1" = "SBT[1]",
  "pub.SBT_2" = "SBT[2]",
  "pub.SBT_3" = "SBT[3]",
  "pub.SBT_4" = "SBT[4]",
  "pub.IEt.CWT" = "IE.CWT",
  "pub.IEt.SBT" = "IE.SBT"
)

trace_data$Coefficient <- factor(trace_data$Coefficient, level = c(unique(trace_data$Coefficient)))
p <- ggplot(trace_data, aes(x = Iter, y = Value, color = type, linewidth = type)) +
  geom_line(alpha = 0.6) +
  facet_grid(Coefficient ~ Method, scales = "free",
             labeller = labeller(Coefficient = coeff_labels,
                                 Method = method_labels,
                                 .default = label_parsed)) +
  # Line sizes and colors
  scale_linewidth_manual(values = c("MCMC" = 0.5, "VB" = 1.5)) +
  scale_color_manual(values = c("MCMC" = "black", "VB" = "black")) +
  # Theme
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    axis.title   = element_text(size = 22),
    axis.text    = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(linewidth = "none", color = "none") +  # Remove legend，linewidth
  labs(
    x = "Iteration",
    y = "Estimated coefficients of JSTVC"
  )

ggsave(p, file = "./figure/FigS20_trace_plots_realData_Tanzania.pdf", width = 18, height = 18)



