source("./LoadPackages/RDependPackages.R")
load("./data/Google_Kenya_Tanzania_Map.RData")
load("./data/Kenya_Score_Data_r.RData")
library(RColorBrewer)

da <- Kenya_Score_Data[Kenya_Score_Data$Year  == 2011, c(1, 3, 6:8, 23)]
p1 <- ggmap(ken.map, darken = c(0, "white")) +
  geom_point(data = da, aes(x = Longitude,
                            y = Latitude,
                            group = as.factor(flag),
                            shape = as.factor(flag),
                            col = (Prevalence)),
             size = 5) +
  scale_shape_manual("Subregions in Kenya ", values = c(20, 18, 17)
                     , label = c("Northwest",
                                 "Northeast",
                                 "South")
  ) +
  scale_alpha_manual("", values = c(0.5, 0.5, 0.5)) +
  coord_fixed(ylim = c(-0.6, 0), xlim = c(34, 35)) +
  scale_x_continuous(limits = c(34, 35),
                     breaks = seq(34, 35, 0.5),
                     labels = paste0(seq(34, 35, 0.5), "° E")) +
  scale_y_continuous(limits = c(-0.6, 0),
                     breaks = seq(-0.6, 0, 0.2),
                     labels = ifelse(seq(-0.6, 0, 0.2) < 0,
                                     paste0(seq(-0.6, 0, 0.2), "° S"),
                                     paste0(seq(-0.6, 0, 0.2), "° N"))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text = element_text(size = 30, colour = "black")
        , axis.title   = element_text(size = 35, colour = "black")
        , legend.title = element_text(size = 35, colour = "black")
        , legend.text  = element_text(size = 35, colour = "black")
        , strip.text   = element_text(size = 35, colour = "black")
        , strip.background = element_rect(colour = "grey100", fill = "grey100")
        , legend.background = element_rect(colour = 'transparent', fill = 'transparent')
        , legend.key.width = unit(3,"line")
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.position  =  c("top")
        , legend.key = element_blank()
        , legend.justification = "right"
  )

int <- 0.1; m <- c(0, 1)
myPalette <- rev(heat.colors(20))
sc <- scale_colour_gradientn(colours = myPalette
                             , limits = m
                             , name = "Relative reduction   "
                             , breaks = c(m[1], round(c(seq(m[1], m[2] - int, int)), 0)[-1], m[2])
                             , labels = c(m[1], round(c(seq(m[1], m[2] - int, int)), 0)[-1], m[2]))

p1 <- p1 + sc + guides(shape = guide_legend(nrow = 1, byrow = T, order = 2, title.position = "top"),
       color = "none")

# p1
ggsave(p1, file = paste0("./figure/Fig2_2_Kenya.pdf"),
       width = 16, height = 8)



load("./data/Tanzania_Score_Data_r.RData")
da <- Tanzania_Score_Data[Tanzania_Score_Data$Year  == 2011, c(1, 3, 6:8, 23)]


p2 <- ggmap(tan.map, darken = c(0, "white")) +
  geom_point(data = da, aes(x = Longitude,
                            y = Latitude,
                            group = as.factor(flag),
                            shape = as.factor(flag),
                            color  = Prevalence),
             size = 5) +
  scale_shape_manual("\n Subregions in Tanzania ", values = c(20, 18, 17)
                     , label = c("West", "East")) +
  coord_fixed(xlim = c(31.5, 34),
              ylim = c(-3, -2)) +
  scale_x_continuous(limits = c(31.5, 34),
                     breaks = seq(31.5, 34, 0.5),
                     labels = paste0(seq(31.5, 34, 0.5), "° E")) +
  scale_y_continuous(limits = c(-3, -2),
                     breaks = seq(-3, -2, 0.5),
                     labels = ifelse(seq(-3, -2, 0.5) < 0,
                                     paste0(seq(-3, -2, 0.5), "° S"),
                                     paste0(seq(-3, -2, 0.5), "° N"))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text = element_text(size = 30, colour = "black")
        , axis.title   = element_text(size = 35, colour = "black")
        , legend.title = element_text(size = 35, colour = "black")
        , legend.text  = element_text(size = 35, colour = "black")
        , strip.text   = element_text(size = 35, colour = "black")
        , strip.background = element_rect(colour = "grey100", fill = "grey100")
        , legend.background = element_rect(colour = 'transparent', fill = 'transparent')
        , legend.key.width = unit(5,"line")
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.position  =  c("top")
        , legend.key = element_blank()
        , legend.spacing.x = unit(100,"pt")
        , legend.justification = "right"
  )

int <- 0.2; m <- c(0, 1)
myPalette <- rev(heat.colors(20))
sc <- scale_colour_gradientn(colours = myPalette
                             , limits = m
                             , name = "Prevalence   "
                             , breaks = c(m[1], round(c(seq(m[1], m[2] - int, int)), 1)[-1], m[2])
                             , labels = c(m[1], round(c(seq(m[1], m[2] - int, int)), 1)[-1], m[2]))

p2 <- p2 + sc + guides(shape = guide_legend(nrow = 1, byrow = T, order = 2, title.position = "top"),
                       color = guide_colorbar(order = 1,title.position = "top" ))
ggsave(p2, file = paste0("./figure/Fig2_1_Tanzania.pdf"),
       width = 16, height = 8)
