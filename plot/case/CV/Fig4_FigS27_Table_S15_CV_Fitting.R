rm(list = ls())
source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
load("./data/Tanzania_Score_Data_r.RData")
Score_Data <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
Score_Data$Prevalence <- log(-log(1 - Score_Data$Prevalence))
setDF(Score_Data)
load(paste0("./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
#-----------------------------------------
set.seed(1234)
#-----------------------------------------
JSTVC <- NULL
for(l in c(1:5))
{
  temp <- apply(CV_Ranalysis$fitted.Pred.ens[[l]], c(1, 2), median) %>% as.data.frame()
  # library(tidyr)
  # library(dplyr)

  temp_long <- temp %>%
    mutate(Year = row_number() + 2010) %>%
    pivot_longer(
      cols = -Year,
      names_to = "Village_ID",
      values_to = "JSTVC"
    ) %>% as.data.frame()
  JSTVC <- rbind(JSTVC, temp_long)
}


#------JSTVC without random effects
X_vars <- c("Intercept"
            , paste0("CWT_", 1:4)
            , paste0("SBT_", 1:4)
            , "IEt.CWT"
            , "IEt.SBT")


str.1 <- paste0("Prevalence ~ -1 + ", paste0(X_vars, collapse = " + "))

str.2 <- paste0("Prevalence ~ -1 + f(Year, model = 'ar1') + f(flag, model = 'iid') + ",
                paste0(X_vars, collapse = " + "))

str.3 <- paste0("Prevalence ~ -1 + f(Year, model = 'ar1') + f(Study_Arm, model = 'iid') + ",
                paste0(X_vars, collapse = " + "))

Null_model <- inla(Prevalence ~ 1,
                   family  = "gaussian",
                   data    = Score_Data,
                   verbose = FALSE,
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE))

JSTVC_xi  <- inla(as.formula(str.1),
                  family  = "gaussian",
                  data    = Score_Data,
                  verbose = FALSE,
                  control.predictor = list(compute = TRUE),
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))
Sub_AR1  <- inla(as.formula(str.2),
                 family  = "gaussian",
                 data    = Score_Data,
                 verbose = FALSE,
                 control.predictor = list(compute = TRUE),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))
Arm_AR1  <- inla(as.formula(str.3),
                 family  = "gaussian",
                 data    = Score_Data,
                 verbose = FALSE,
                 control.predictor = list(compute = TRUE),
                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))

Score_Data$Fitted.JSTVC_xi <- (1 - exp(-exp(JSTVC_xi[["summary.fitted.values"]][["0.5quant"]])))*1
Score_Data$Fitted.Sub_AR1  <- (1 - exp(-exp(Sub_AR1[["summary.fitted.values"]][["0.5quant"]])))*1
Score_Data$Fitted.Arm_AR1  <- (1 - exp(-exp(Arm_AR1[["summary.fitted.values"]][["0.5quant"]])))*1
Score_Data$Fitted.Null     <- (1 - exp(-exp(Null_model[["summary.fitted.values"]][["0.5quant"]])))*1

Score_Data$Prevalence <- (1 - exp(- exp(Score_Data$Prevalence)))*1
Fitted <- JSTVC %>% left_join(Score_Data[, c("Village_ID",
                                             "Year",
                                             "Prevalence",
                                             "Fitted.Null",
                                             "Fitted.JSTVC_xi",
                                             "Fitted.Sub_AR1",
                                             "Fitted.Arm_AR1")],
                              by = c("Village_ID", "Year"))

# In-sample goodness-of-fit benchmark against an intercept-only null model
null.mean <- Fitted$Fitted.Null
null.sse <- sum((Fitted$Prevalence - Fitted$Fitted.Null)^2, na.rm = TRUE)

calc_fit_metrics <- function(obs, pred, model.name) {
  sse.model <- sum((obs - pred)^2, na.rm = TRUE)
  data.frame(
    Model = model.name,
    RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    MAE = mean(abs(obs - pred), na.rm = TRUE),
    Fitted_R2 = 1 - sse.model / null.sse,
    DIC = NA_real_,
    WAIC = NA_real_
  )
}

Fitted.GOF.Table <- rbind(
  data.frame(
    Model = "Null",
    RMSE = sqrt(mean((Fitted$Prevalence - Fitted$Fitted.Null)^2, na.rm = TRUE)),
    MAE = mean(abs(Fitted$Prevalence - Fitted$Fitted.Null), na.rm = TRUE),
    Fitted_R2 = 0,
    DIC = Null_model$dic$dic,
    WAIC = Null_model$waic$waic
  ),
  calc_fit_metrics(Fitted$Prevalence, Fitted$JSTVC, "JSTVC"),
  calc_fit_metrics(Fitted$Prevalence, Fitted$Fitted.JSTVC_xi, "JSTVC_xi"),
  calc_fit_metrics(Fitted$Prevalence, Fitted$Fitted.Arm_AR1, "Arm_AR1"),
  calc_fit_metrics(Fitted$Prevalence, Fitted$Fitted.Sub_AR1, "Sub_AR1")
) %>%
  mutate(
    DIC = c(
      Null_model$dic$dic,
      CV_Ranalysis$Goodness.of.fit$Overall$DIC,
      JSTVC_xi$dic$dic,
      Arm_AR1$dic$dic,
      Sub_AR1$dic$dic
    ),
    WAIC = c(
      Null_model$waic$waic,
      CV_Ranalysis$Goodness.of.fit$Overall$WAIC,
      JSTVC_xi$waic$waic,
      Arm_AR1$waic$waic,
      Sub_AR1$waic$waic
    )
  ) %>%
  mutate(across(c(RMSE, MAE, Fitted_R2, DIC, WAIC), ~ round(.x, 4)))

Fitted.GOF.Table <- Fitted.GOF.Table[match(c("Null", "JSTVC_xi", "Arm_AR1", "Sub_AR1", "JSTVC"),
                                           Fitted.GOF.Table$Model), ]

Fitted.GOF.Table$RMSE <- NULL
Fitted.GOF.Table$MAE <- NULL
writexl::write_xlsx(
  list(Fitted_GOF_Benchmark = as.data.frame(Fitted.GOF.Table)),
  path = "./result/summary/Table_S15_good_of_fitting.xlsx"
)

Fitted.GOF.Table



load(file = "./result/case/SCORE_CV_models_RMSE_CRPS_ie.RData")
Res <- cv.Pred.details[!is.na(cv.Pred.details$true.Prevalence),]


pdf(file  = paste0("./figure/Fig4_Obser_Pred.pdf"), width = 16, height = 8)
par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.8, 0))
par(mfrow = c(1, 2),
    cex      = 1.5,
    cex.axis = 1.2,
    cex.lab  = 1.5,
    cex.main = 1.5,
    lwd      = 2)
plot(Res$true.Prevalence*1,
     Res$pred.JSTVC*1,
     xlim =  c(0, 1),
     ylim =  c(0, 1),
     xlab = "Observed prevalence",
     ylab = "Predicted prevalence",
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5))
points(Res$true.Prevalence*1, Res$pred.JSTVC_xi*1, col= adjustcolor("blue", alpha.f = 1), pch = 1)
points(Res$true.Prevalence*1, Res$pred.Arm_AR1*1, col= adjustcolor("green", alpha.f = 0.5), pch = 8)
points(Res$true.Prevalence*1, Res$pred.Sub_AR1*1, col= adjustcolor("red", alpha.f = 0.5), pch = 3)

legend("bottomright",
       legend = c(expression(JSTVC[-xi]),
                  TeX("$Arm-AR_1$"),
                  TeX("$Sub-AR_1$"),
                  "JSTVC"),
       col    = c(adjustcolor("blue",  alpha.f = 1),
                  adjustcolor("red",   alpha.f = 0.5),
                  adjustcolor("green", alpha.f = 0.5),
                  adjustcolor("black", alpha.f = 0.5)),
       inset = c(0, -0.02),
       pch = c(1, 8, 3, 19),
       bty = "n",
       cex = 0.9)
mtext("(A)", side = 3, line = 0.85, adj = -0.1, cex = 2, font = 2)
plot(Fitted$Prevalence*1,
     Fitted$JSTVC*1,
     xlim =  c(0, 1),
     ylim =  c(0, 1),
     xlab = "Observed prevalence",
     ylab = "Fitted prevalence",
     pch  = 19,
     col  = adjustcolor("black", alpha.f = 0.5))
points(Fitted$Prevalence*1, Fitted$Fitted.JSTVC_xi*1, col= adjustcolor("blue", alpha.f = 1), pch = 1)
points(Fitted$Prevalence*1, Fitted$Fitted.Sub_AR1*1, col= adjustcolor("green", alpha.f = 0.5), pch = 8)
points(Fitted$Prevalence*1, Fitted$Fitted.Arm_AR1*1, col= adjustcolor("red", alpha.f = 0.5), pch = 3)

legend("bottomright",
       legend = c(expression(JSTVC[-xi]),
                  TeX("$Arm-AR_1$"),
                  TeX("$Sub-AR_1$"),
                  "JSTVC"),
       col    = c(adjustcolor("blue", alpha.f = 1),
                  adjustcolor("red", alpha.f = 0.5),
                  adjustcolor("green", alpha.f = 0.5),
                  adjustcolor("black", alpha.f = 0.5)),
       inset = c(0, -0.02),
       pch = c(1, 8, 3, 19),
       bty = "n",
       cex = 1)
mtext("(B)", side = 3, line = 0.85, adj = -0.1, cex = 2, font = 2)
dev.off()


# Spatiotemporal residual maps for in-sample fitted values
Residual.map <- Fitted %>%
  left_join(
    unique(Score_Data[, c("Village_ID", "Year", "Longitude", "Latitude", "flag")]),
    by = c("Village_ID", "Year")
  ) %>%
  transmute(
    Village_ID,
    YEAR = Year,
    LON = Longitude,
    LAT = Latitude,
    flag,
    JSTVC = Prevalence - JSTVC,
    `JSTVC-xi` = Prevalence - Fitted.JSTVC_xi,
    `Sub-AR1` = Prevalence - Fitted.Sub_AR1,
    `Arm-AR1` = Prevalence - Fitted.Arm_AR1
  ) %>%
  pivot_longer(
    cols = c("JSTVC", "JSTVC-xi", "Sub-AR1", "Arm-AR1"),
    names_to = "Model",
    values_to = "Residual"
  ) %>%
  mutate(
    Country = ifelse(grepl("^KEN", Village_ID), "Kenya", "Tanzania"),
    Model = factor(Model, levels = c("JSTVC-xi", "Arm-AR1", "Sub-AR1", "JSTVC")),
    YEAR = factor(YEAR, levels = 2011:2015)
  )

Residual.range <- max(abs(Residual.map$Residual), na.rm = TRUE)


#Comparisons of outliers

temp <- as.data.frame(Residual.map[Residual.map$Village_ID %in% c("KEN212", "KEN199", "KEN086"), ])
aggregate(Residual ~ Village_ID,
          data = temp,
          FUN = mean)




load("./data/Google_Kenya_Tanzania_Map.RData")

model_labeller <- c(
  "JSTVC-xi" = "\"JSTVC\"[xi]",
  "Arm-AR1" = "\"Arm-AR\"[1]",
  "Sub-AR1" = "\"Sub-AR\"[1]",
  "JSTVC" = "\"JSTVC\""
)

plot_residual_map <- function(map.data,
                              plot.data,
                              shape.values,
                              shape.labels,
                              shape.title,
                              xlim,
                              ylim,
                              xbreaks,
                              ybreaks) {
  ggmap(map.data, darken = c(0, "white")) +
    geom_point(data = plot.data,
               aes(x = LON,
                   y = LAT,
                   shape = flag,
                   colour = Residual),
               size = 1.7) +
    facet_grid(
      Model ~ YEAR,
      space = "free",
      labeller = labeller(Model = as_labeller(model_labeller, label_parsed))
    ) +
    coord_fixed(xlim = xlim, ylim = ylim) +
    scale_shape_manual(shape.title,
                       values = shape.values,
                       labels = shape.labels) +
    scale_colour_gradient2(
      low = "#2b6cb0",
      mid = "white",
      high = "#c53030",
      midpoint = 0,
      limits = range(na.omit(Residual.tanzania$Residual)),
      name = "Residuals (observed - predicted)"
    ) +
    scale_x_continuous(
      limits = xlim,
      breaks = xbreaks,
      labels = paste0(xbreaks, "° E")
    ) +
    scale_y_continuous(
      limits = ylim,
      breaks = ybreaks,
      labels = ifelse(
        ybreaks < 0,
        paste0(abs(ybreaks), "° S"),
        paste0(ybreaks, "° N")
      )
    ) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(
      axis.text = element_text(size = 16, colour = "black"),
      axis.title = element_text(size = 18, colour = "black"),
      legend.title = element_text(size = 18, colour = "black"),
      legend.text = element_text(size = 16, colour = "black"),
      strip.text = element_text(size = 16, colour = "black"),
      strip.background = element_rect(colour = "grey100", fill = "grey100"),
      legend.background = element_rect(colour = "transparent", fill = "transparent"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.key = element_blank(),
      legend.key.width = unit(3.8, "line")
    ) +
    guides(
      shape = guide_legend(
        nrow = 1,
        byrow = TRUE,
        order = 1,
        title.position = "top",
        keywidth  = unit(2.8, "lines"),
        override.aes = list(size = 5)
      ),
      colour = guide_colorbar(order = 2, title.position = "top")
    )
}
Residual.map <- Residual.map[!is.na(Residual.map$Residual),]
Residual.kenya <- Residual.map[Residual.map$Country == "Kenya", ]
Residual.kenya$flag <- factor(
  Residual.kenya$flag,
  levels = c("Northwestern", "Northeastern", "Southern")
)

Residual.tanzania <- Residual.map[Residual.map$Country == "Tanzania", ]
Residual.tanzania$flag <- factor(
  Residual.tanzania$flag,
  levels = c("Western", "Eastern")
)

p.kenya <- plot_residual_map(
  map.data = ken.map,
  plot.data = Residual.kenya,
  shape.values = c(20, 18, 17),
  shape.labels = c("Northwest", "Northeast", "South"),
  shape.title = "Subregions in Kenya",
  xlim = c(34, 35),
  ylim = c(-0.6, 0),
  xbreaks = seq(34, 35, 0.4),
  ybreaks = seq(-0.6, 0, 0.2)
)

p.tanzania <- plot_residual_map(
  map.data = tan.map,
  plot.data = Residual.tanzania,
  shape.values = c(20, 17),
  shape.labels = c("West", "East"),
  shape.title = "Subregions in Tanzania",
  xlim = c(31.5, 34),
  ylim = c(-3, -2),
  xbreaks = seq(31.8, 33.8, 1),
  ybreaks = seq(-3, -2, 0.5)
)

pdf("./figure/FigS27_1_kenya.pdf", width = 18, height = 10)
print(p.kenya)
dev.off()

pdf("./figure/FigS27_2_tanzania.pdf", width = 18, height = 10)
print(p.tanzania)
dev.off()







