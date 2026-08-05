rm(list = ls())
source("./LoadPackages/RDependPackages.R")
load(paste0("./result/single_smoothed/sim_smoothed_MCMC_FALSE_250.RData"))
elapsed.vb.1 <- round(abs(run.time)[3]/3600, 1)
JSTVC.vb   <- CV_Ranalysis[[1]]$Process.monitoring$VB.EnKS.Iter#[, -2]
vb_long    <- JSTVC.vb %>%
              pivot_longer(
                cols = starts_with("pub."),
                names_to = "Coefficient",
                values_to = "Value"
              ) %>%
              mutate(Method = "VB (n = 250)")



load(paste0("./result/single_smoothed/sim_smoothed_MCMC_TRUE_250.RData"))
elapsed.mcmc.1 <- round(abs(run.time)[3]/3600, 1)
JSTVC.mcmc <- CV_Ranalysis[[1]]$Process.monitoring$VB.EnKS.Iter#[, -2]
mcmc_long  <- JSTVC.mcmc %>%
  filter(Iter >= 0e4) %>%
  pivot_longer(
    cols = starts_with("pub."),
    names_to = "Coefficient",
    values_to = "Value"
  ) %>%
  mutate(Method = "MCMC (n = 250)")

vb_long$type   <- "VB"
mcmc_long$type <- "MCMC"
trace_data.250 <- bind_rows(vb_long, mcmc_long)
#------------------------------n = 298----------------------------------------
load(paste0("./result/single_smoothed/sim_smoothed_MCMC_FALSE_298.RData"))
elapsed.vb.2 <- round(abs(run.time)[3]/3600, 1)
JSTVC.vb <- CV_Ranalysis$Process.monitoring$VB.EnKS.Iter#[, -2]
vb_long  <- JSTVC.vb %>%
  pivot_longer(
    cols = starts_with("pub."),
    names_to = "Coefficient",
    values_to = "Value"
  ) %>%
  mutate(Method = "VB (n = 298)")



load(paste0("./result/single_smoothed/sim_smoothed_MCMC_TRUE_298.RData"))
elapsed.mcmc.2 <- round(abs(run.time)[3]/3600, 1)
JSTVC.mcmc <- CV_Ranalysis$Process.monitoring$VB.EnKS.Iter#[, -2]
mcmc_long  <- JSTVC.mcmc %>%
  filter(Iter >= 0e4) %>%
  pivot_longer(
    cols       = starts_with("pub."),
    names_to   = "Coefficient",
    values_to  = "Value"
  ) %>%
  mutate(Method = "MCMC (n = 298)")

vb_long$type   <- "VB"
mcmc_long$type <- "MCMC"
trace_data <- bind_rows(trace_data.250, vb_long, mcmc_long)
############################################################
############################################################
method_labels <- c(
  "MCMC (n = 250)" = paste0("MCMC~'(n = 250; running time = ", elapsed.mcmc.1, "h)'"),
  "VB (n = 250)"   = paste0("VB~'(n = 250; running time = ", elapsed.vb.1, "h)'"),
  "MCMC (n = 298)" = paste0("MCMC~'(n = 298; running time = ", elapsed.mcmc.2, "h)'"),
  "VB (n = 298)"   = paste0("VB~'(n = 298; running time = ", elapsed.vb.2, "h)'")
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
trace_data$Method <- factor(trace_data$Method, level = c("MCMC (n = 250)", "VB (n = 250)",
                                                         "MCMC (n = 298)", "VB (n = 298)"))

p <- ggplot(trace_data, aes(x = Iter, y = Value, color = type, linewidth = type)) +
  geom_line(alpha = 0.6) +
  geom_hline(data = subset(trace_data, Coefficient == "pub.Intercept"),
             aes(yintercept = 5), linetype = "dashed", color = "blue") +
  geom_hline(data = subset(trace_data, Coefficient != "pub.Intercept"),
             aes(yintercept = -1), linetype = "dashed", color = "blue") +
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
    strip.text.y = element_text(size = 22),
    axis.title   = element_text(size = 25),
    axis.text    = element_text(size = 20),
    axis.text.y    = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(linewidth = "none", color = "none") +  # Remove legend，linewidth
  labs(
    x = "Iteration",
    y = "Trace of estimated coefficients of JSTVC"
  )

ggsave(p, file = "./figure/FigS8_smoothed_trace_plots.pdf", width = 22, height = 18)



