rm(list = ls())
#-----------------------------------------
source("./LoadPackages/RDependPackages.R")
load("./data/Google_Kenya_Tanzania_Map.RData")
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site
load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site

#-----------------------------------------
p1 <- ggmap(tan.map, darken = c(0, "white")) +
  geom_point(data = Tan.Site, aes(x = LON,
                              y = LAT,
                              group = as.factor(Study_Arm),
                              col = as.factor(Study_Arm)),
             size = 3) +
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
  labs(color = "Arms") +
  theme(axis.text = element_text(size = 23, colour = "black")
        , axis.title   = element_text(size = 28, colour = "black")
        , legend.title = element_text(size = 25, colour = "black")
        , legend.text  = element_text(size = 25, colour = "black")
        , strip.text   = element_text(size = 25, colour = "black")
        , strip.background = element_rect(colour = "grey100", fill = "grey100")
        , legend.background = element_rect(colour = 'transparent', fill = 'transparent')
        , legend.key.width = unit(5,"line")
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.position  =  c("top")
  )  +
  guides(col = guide_legend(override.aes = list(size = 5),
                            nrow = 1, byrow = TRUE))

p2 <- ggmap(ken.map, darken = c(0, "white")) +
  geom_point(data = Ken.Site, aes(x = LON,
                              y = LAT,
                              group = as.factor(Study_Arm),
                              col = as.factor(Study_Arm)),
             size = 3) +
  coord_fixed(ylim = c(-0.6, 0), xlim = c(34, 35)) +
  scale_x_continuous(limits = c(34, 35),
                     breaks = seq(34, 35, 0.2),
                     labels = paste0(seq(34, 35, 0.2), "° E")) +
  scale_y_continuous(limits = c(-0.6, 0),
                     breaks = seq(-0.6, 0, 0.2),
                     labels = ifelse(seq(-0.6, 0, 0.2) < 0,
                                     paste0(seq(-0.6, 0, 0.2), "° S"),
                                     paste0(seq(-0.6, 0, 0.2), "° N"))) +
  xlab("Longitude") + ylab("Latitude") +
  labs(color = "Arms") +
  theme(axis.text = element_text(size = 23, colour = "black")
        , axis.title   = element_text(size = 28, colour = "black")
        , legend.title = element_text(size = 25, colour = "black")
        , legend.text  = element_text(size = 25, colour = "black")
        , strip.text   = element_text(size = 25, colour = "black")
        , strip.background = element_rect(colour = "grey100", fill = "grey100")
        , legend.background = element_rect(colour = 'transparent', fill = 'transparent')
        , legend.key.width = unit(5,"line")
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.position  =  c("top")
  )  +
  guides(col = guide_legend(override.aes = list(size = 5), nrow = 1, byrow = TRUE))


ggsave(p1, file = paste0("./figure/FigS1_1_Tanania_Arms.pdf"),
      width = 15, height = 8)

ggsave(p2, file = paste0("./figure/FigS1_2_Kenya_Arms.pdf"),
       width = 15, height = 8)
