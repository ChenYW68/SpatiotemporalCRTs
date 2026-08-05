rm(list = ls())
source("./LoadPackages/RDependPackages.R")
load(paste0("./result/single_smoothed/sim_smoothed_MCMC_TRUE_298.RData"))
true_df <- setDF(sim_Data)  %>%
  dplyr::select(time.index, LON, LAT, W_ts) %>%
  rename(Value = W_ts) %>%
  mutate(Method = "Simulated") %>% setorderv(c("LON", "LAT", "time.index"))

mcmc_df  <- setDF(CV_Ranalysis$Process.monitoring$GRFs.process) %>%
  filter(iter > 1.5e4) %>%
  dplyr::group_by(time.index, LON, LAT) %>%
  dplyr::summarise(Value = mean(intercept.GRF), .groups = "drop") %>%
  as.data.frame() %>%
  mutate(Method = "MCMC")%>% setorderv(c("LON", "LAT", "time.index"))


load(paste0("./result/single_smoothed/sim_smoothed_MCMC_FALSE_298.RData"))
vb_df <- setDF(CV_Ranalysis$Process.monitoring$GRFs.process) %>%
  dplyr::select(time.index, LON, LAT, intercept.GRF) %>%
  rename(Value = intercept.GRF) %>%
  mutate(Method = "VB") %>% setorderv(c("LON", "LAT", "time.index"))

head(true_df)
head(mcmc_df)
head(vb_df)


plot_data <- bind_rows(true_df, mcmc_df, vb_df)

# Make Method an ordered factor so rows are True -> MCMC -> VB
plot_data$Method <- factor(plot_data$Method, levels = c("Simulated", "MCMC", "VB"))

# Plot
# Define color limits across all data for consistent scale
val_range <- c(floor(min(plot_data$Value)*10)/10, ceiling(max(plot_data$Value)*10)/10)
val_mid   <- mean(plot_data$Value, na.rm = TRUE)

# Plot
time_labels <- c(
  "1" = "time = 1",
  "2" = "time = 2",
  "3" = "time = 3",
  "4" = "time = 4",
  "5" = "time = 5")

p <- ggplot(plot_data[plot_data$LAT > -1, ], aes(x = LON, y = LAT, col = Value)) +
  geom_tile() +
  geom_point(size = 1.2) +
  facet_grid(vars(Method), vars(time.index),
             labeller = labeller(time.index = time_labels)) +

  scale_x_continuous(limits = c(34, 35),
                     breaks = c(34, 34.4, 34.8),
                     labels = paste0(c(34, 34.4, 34.8), "° E")) +
  scale_y_continuous(limits = c(-0.6, 0),
                     breaks = seq(-0.6, 0, 0.2),
                     labels = ifelse(seq(-0.6, 0, 0.2) < 0,
                                     paste0(seq(-0.6, 0, 0.2), "° S"),
                                     paste0(seq(-0.6, 0, 0.2), "° N"))) +

  theme_light() +
  theme(
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 18, colour = "black"),
    legend.title = element_text(size = 18, colour = "black"),
    legend.text  = element_text(size = 14, colour = "black"),
    strip.text   = element_text(size = 16, colour = "black"),
    strip.background = element_rect(colour = "grey80", fill = "grey80"),
    legend.background = element_rect(colour = 'transparent', fill = 'transparent'),
    legend.key.width = unit(10,"line"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  ) +
  labs(
    x = 'Longitude',
    y = 'Latitude',
    col = TeX("Surface: $xi_t(s)$")
  )+
  guides(
    color = guide_colorbar(
      title.position = "left",  # moves title above the bar
      title.vjust = 1        # centers the title over the bar
    )
  )

int <- 0.5
library(RColorBrewer)
#https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html
# display.brewer.all(colorblindFriendly=TRUE)
m <- c(-1, 1) #rev(heat.colors(20))#
myPalette <- colorRampPalette(rev(brewer.pal(20, "Spectral"))) #"Spectral" #RdYlGn
sc <- scale_colour_gradientn(colours = myPalette(100) #myPalette#
                             , limits = m
                             # , name = "log(Prevalence) and estimated effects   "
                             , breaks = c(seq(m[1], m[2], by = int))
                             , labels = c(seq(m[1], m[2], by = int)))

p <- p + sc

ggsave(p, file = paste0("./figure/FigS10_smoothed_surface.pdf"),
       width = 14, height = 8)
