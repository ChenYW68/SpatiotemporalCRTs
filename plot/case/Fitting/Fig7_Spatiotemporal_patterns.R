# beta <- data.frame(Variable = rownames(beta.1$mu.beta),
#                    beta = round(as.vector(beta.1$mu.beta[, 1]), 2),
#                    sd = round(beta.1$cov.prob[, n.c - 2], 2),
#                    # post.prob = beta.1$cov.prob[, n.c],
#                    group = "All regions",
#                    model = "Non-joint models")
#
rm(list = ls())
source("./LoadPackages/RDependPackages.R")

inverse_cloglog <- function(eta) {
  1 - exp(-exp(eta))
}

load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site

load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site


# total.G <- as.matrix(bdiag(Ken.G, G.mat))
S1 <- Ken.Site[Ken.Site$flag =="Northwestern", ]
S2 <- Ken.Site[Ken.Site$flag =="Northeastern", ]
S3 <- Ken.Site[Ken.Site$flag =="Southern", ]
S4 <- Tan.Site[Tan.Site$flag =="Western", ]
S5 <- Tan.Site[Tan.Site$flag =="Eastern", ]
Site <- rbind(S1, S2, S3, S4, S5)
Site$flag <- as.character(Site$flag)
ind.1 <- which(Site$flag =="Northwestern")
ind.2 <- which(Site$flag =="Northeastern")
ind.3 <- which(Site$flag =="Southern")
ind.4 <- which(Site$flag =="Western")
ind.5 <- which(Site$flag =="Eastern")

Score_Data       <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
load(paste0("./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
beta.1 <- CV_Ranalysis$update.Para.List[[1]]$beta

mean.RF.1 <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[1]][,1:length(ind.1),], c(1, 2), quantile, prob = c(0.025, 0.5, 0.975))
sd.RF.1 <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[1]][,ind.1,], 1:2, sd)
range(sd.RF.1)

ind.leng <- (1 + length(ind.1)): (length(ind.1) + length(ind.2))

mean.RF.2 <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[2]][, ind.leng,], c(1, 2), quantile, prob = c(0.025, 0.5, 0.975))
sd.RF.2   <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[2]][,ind.leng,], 1:2, sd)
range(sd.RF.2)


ind.leng  <- (1 + length(ind.1) + length(ind.2)): (length(ind.1) + length(ind.2) + length(ind.3))
mean.RF.3 <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[3]][,ind.leng,], c(1, 2), quantile, prob = c(0.025, 0.5, 0.975))
sd.RF.3   <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[3]][,ind.leng,], 1:2, sd)
range(sd.RF.3)


ind.leng  <- (1 + length(ind.1) + length(ind.2) + length(ind.3)):(length(ind.1) + length(ind.2) + length(ind.3)+ length(ind.4))
mean.RF.4 <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[4]][,ind.leng,], c(1, 2), quantile, prob = c(0.025, 0.5, 0.975))
sd.RF.4   <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[4]][,ind.leng,], 1:2, sd)
range(sd.RF.4)

ind.leng  <- (1 + length(ind.1) + length(ind.2) + length(ind.3) + length(ind.4)): (length(ind.1) + length(ind.2) + length(ind.3) + length(ind.4)+ length(ind.5))
mean.RF.5 <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[5]][,ind.leng,], c(1, 2), quantile, prob = c(0.025, 0.5, 0.975))
sd.RF.5   <- apply(CV_Ranalysis$fitted.slope.Hv.Zg.ens[[5]][,ind.leng,], 1:2, sd)
range(sd.RF.5)




sRF.1     <- as.data.frame(mean.RF.1[2,,])
Col.Name  <- colnames(sRF.1)
sRF.1$YEAR <- 2011:2015
sRF.1 <- sRF.1 %>% pivot_longer(cols = Col.Name,
                                names_to ='Village_ID',
                                values_to ='Wts') %>% as.data.frame()


sRF.2     <- as.data.frame(mean.RF.2[2,,])
Col.Name  <- colnames(sRF.2)
sRF.2$YEAR <- 2011:2015
sRF.2 <- sRF.2 %>% pivot_longer(cols = Col.Name,
                                names_to ='Village_ID',
                                values_to ='Wts') %>% as.data.frame()


sRF.3     <- as.data.frame(mean.RF.3[2,,])
Col.Name  <- colnames(sRF.3)
sRF.3$YEAR <- 2011:2015
sRF.3 <- sRF.3 %>% pivot_longer(cols = Col.Name,
                                names_to ='Village_ID',
                                values_to ='Wts') %>% as.data.frame()


sRF.4     <- as.data.frame(mean.RF.4[2,,])
Col.Name  <- colnames(sRF.4)
sRF.4$YEAR <- 2011:2015
sRF.4 <- sRF.4 %>% pivot_longer(cols = Col.Name,
                                names_to ='Village_ID',
                                values_to ='Wts') %>% as.data.frame()

sRF.5     <- as.data.frame(mean.RF.5[2,,])
Col.Name  <- colnames(sRF.5)
sRF.5$YEAR <- 2011:2015
sRF.5 <- sRF.5 %>% pivot_longer(cols = Col.Name,
                                names_to ='Village_ID',
                                values_to ='Wts') %>% as.data.frame()



sd.RF.1     <- as.data.frame(sd.RF.1)
Col.Name  <- colnames(sd.RF.1)
sd.RF.1$YEAR <- 2011:2015
sd.RF.1 <- sd.RF.1 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='Wts') %>% as.data.frame()


sd.RF.2     <- as.data.frame(sd.RF.2)
Col.Name  <- colnames(sd.RF.2)
sd.RF.2$YEAR <- 2011:2015
sd.RF.2 <- sd.RF.2 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='Wts') %>% as.data.frame()


sd.RF.3     <- as.data.frame(sd.RF.3)
Col.Name  <- colnames(sd.RF.3)
sd.RF.3$YEAR <- 2011:2015
sd.RF.3 <- sd.RF.3 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='Wts') %>% as.data.frame()

sd.RF.4     <- as.data.frame(sd.RF.4)
Col.Name  <- colnames(sd.RF.4)
sd.RF.4$YEAR <- 2011:2015
sd.RF.4 <- sd.RF.4 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='Wts') %>% as.data.frame()


sd.RF.5     <- as.data.frame(sd.RF.5)
Col.Name  <- colnames(sd.RF.5)
sd.RF.5$YEAR <- 2011:2015
sd.RF.5 <- sd.RF.5 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='Wts') %>% as.data.frame()


sRF <- rbind(sRF.1, sRF.2, sRF.3, sRF.4, sRF.5) %>% left_join(Site[, c(1:5)], by = c("Village_ID"))

sRF.Sd <- rbind(sd.RF.1, sd.RF.2, sd.RF.3, sd.RF.4, sd.RF.5) %>% left_join(Site[, c(1:5)], by = c("Village_ID"))
setDF(Score_Data)

Da.prevalence <- Score_Data[, c(2, 1, 3, 5, 6:8)] #"Year", Village_ID, Prevalence, Study_Arm, Longitude, Latitude
Da.treatment.effect <- Score_Data[, c(2, 1, 5, 6:8)]
Da.treatment.effect$Eta.fix <-
  as.matrix(Score_Data[, c(9:16, 19:20)]) %*%
    matrix(as.vector(beta.1$mu.beta[-1, 1]), ncol = 1)
Da.treatment.effect$Wts <- #Da.prevalence$Prevalence -
                            (inverse_cloglog(beta.1$mu.beta[1, 1] + Da.treatment.effect$Eta.fix + sRF$Wts) -
                            inverse_cloglog(beta.1$mu.beta[1, 1] + sRF$Wts))

colnames(Da.prevalence) <- colnames(sRF)[1:7]
setnames(Da.treatment.effect, c("Longitude", "Latitude", "Year"), c("LON", "LAT", "YEAR"))



# range(Da$Wts, na.rm = T)
Da.prevalence$Group       <- "Prevalence"
Da.treatment.effect$Group <- "Treatment effects"
sRF$Group                 <- "Random effects"

sRF <- sRF %>% left_join(Da.treatment.effect[, c("Village_ID", "YEAR", "Eta.fix")], by = c("Village_ID", "YEAR"))




sRF$Wts <-
  inverse_cloglog(beta.1$mu.beta[1, 1] + sRF$Wts)
sRF <- sRF[, 1:8]
Da.prevalence <- Da.prevalence[, c("YEAR", "Village_ID", "Wts", "Study_Arm", "LON", "LAT", "flag", "Group")]
Da.treatment.effect <- Da.treatment.effect[, c("YEAR", "Village_ID", "Wts", "Study_Arm", "LON", "LAT", "flag", "Group")]
Da.sRF <- rbind(sRF, Da.prevalence, Da.treatment.effect)

Da.sRF$Group <- ordered(Da.sRF$Group, level = c("Prevalence", "Treatment effects", "Random effects"))

Da.sRF$flag <- ordered(Da.sRF$flag, level = c("Northwestern", "Northeastern", "Southern", "Western", "Eastern"))
sRF.Sd$flag <- ordered(sRF.Sd$flag, level = c("Northwestern", "Northeastern", "Southern", "Western", "Eastern"))


#Kenya ----
Wts.plot <- Da.sRF[Da.sRF$flag %in% c("Northwestern", "Northeastern", "Southern"), ]

range(sRF$Wts, na.rm = T)
range(Wts.plot$Wts, na.rm = T)
Wts.plot <- Wts.plot[!is.na(Wts.plot$Wts), ]

load("./data/Google_Kenya_Tanzania_Map.RData")
p11 <- suppressMessages({ggmap(ken.map, darken = c(0, "white")) +
  geom_point(data = Wts.plot, aes(x = LON,
                                  y = LAT,
                                  group = as.factor(flag),
                                  shape = as.factor(flag),
                                  col = Wts),
             size = 2) +
  scale_shape_manual("Subregions in Kenya "
                     , values = c(20, 18, 17)
                     , label = c("Northwest",
                                 "Northeast",
                                 "South")
  ) +
  facet_grid(Group ~ YEAR, space = "free") +
  scale_alpha_manual("", values = c(0.5, 0.5, 0.5)) +
  coord_fixed(ylim = c(-0.6, 0), xlim = c(34, 35)) +
  scale_x_continuous(limits = c(34, 35),
                     breaks = seq(34, 35, 0.4),
                     labels = paste0(seq(34, 35, 0.4), "° E")) +
  scale_y_continuous(limits = c(-0.6, 0),
                     breaks = seq(-0.6, 0, 0.2),
                     labels = ifelse(seq(-0.6, 0, 0.2) < 0,
                                     paste0(seq(-0.6, 0, 0.2), "° S"),
                                     paste0(seq(-0.6, 0, 0.2), "° N"))) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text = element_text(size = 16, colour = "black")
        # ,axis.text.x = element_text(hjust = 0.25, size = 35, colour = "black")
        , axis.title   = element_text(size = 18, colour = "black")
        , legend.title = element_text(size = 18, colour = "black")
        , legend.text  = element_text(size = 18, colour = "black")
        , strip.text   = element_text(size = 18, colour = "black")
        # , legend.title = element_blank()
        , strip.background = element_rect(colour = "grey100", fill = "grey100")
        , legend.background = element_rect(colour = 'transparent', fill = 'transparent')
        , legend.key.width = unit(8,"line")
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.position  =  c("top")
        , legend.key = element_blank()
        , legend.spacing.x = unit(50,"pt")
        , legend.justification = "center"
        # , legend.spacing.y = unit(-10.5, "cm")
        , legend.box.spacing = unit(-0.1, "cm")
        # , legend.margin = margin(t = -0.1, unit='cm')
  )})
int <- 1
library(RColorBrewer)
#https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html
# display.brewer.all(colorblindFriendly=TRUE)
m <- c(-0.8, 1)#range(Wts.plot$Wts, na.rm = TRUE)
legend_breaks <- pretty(m, n = 8)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral"))) #"Spectral" #RdYlGn
sc <- scale_colour_gradientn(colours = myPalette(nrow(sRF)) #myPalette#
                             , limits = m
                             , name = "Prevalence and estimated effects"
                             , breaks = legend_breaks
                             , labels = legend_breaks)

library(viridisLite)

sc <- scale_colour_gradientn(
  colours = turbo(256),
  limits = m,
  name = "Prevalence and estimated effects",
  breaks = legend_breaks,
  labels = legend_breaks
)

p11 <- p11 + sc + guides(shape = guide_legend(nrow = 1, byrow = T, order = 2,
                                              title.position = "top" ,
                                              override.aes = list(size = 5),
                                              keywidth  = unit(0.9, "lines")#,
                                              # keyheight = unit(0.8, "lines"),
                                              # legend.background = element_rect(colour = 'transparent', fill = 'transparent')
                                              ),
                         color = guide_colorbar(order = 1, title.position = "top" ))
# p11
ggsave(p11, file = paste0("./figure/Fig7_Kenya_Wts.pdf"),
       width = 18, height = 8)



