# install.packages("./LoadPackages//HDCM_1.0.zip", repos = NULL, type = "win.binary")
rm(list = ls())
#-----------------------------------------
source("./LoadPackages/RDependPackages.R")
source(normalizePath("./JSTVC/R/regCreateGridm.R"))
source(normalizePath("./JSTVC/R/Partitioning.Datasets.R"))
source(normalizePath("./JSTVC/R/Construct.Fixed.effect.Data.R"))
source(normalizePath("./JSTVC/R/util.R"))
source(normalizePath("./JSTVC/R/VB.R"))
source(normalizePath("./JSTVC/R/VB_EnKF.R"))
source(normalizePath("./JSTVC/R/EnKF.R"))
source(normalizePath("./JSTVC/R/JSTVCs.R"))
#-----------------------------------------
Tab <- paste0("./result/case/")
if (!dir.exists(Tab)) {
  dir.create(Tab, recursive = TRUE)
}
#-----------------------------------------
#1. Standard Setting of JSTVC
#-----------------------------------------
X_vars   <- c("Intercept", paste0("CWT_", 1:4), paste0("SBT_", 1:4), "IEt.CWT", "IEt.SBT")
Tab.Name <- paste0(Tab, "/SCORE_CV_models_RMSE_CRPS_ie.RData")
#-----------------------------------------
# Load data
#-----------------------------------------
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site
Ken.G    <- G.mat
load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Tan.G    <- G.mat
Site     <- rbind(Ken.Site, Tan.Site)
#-----------------------------------------
# Define indices for regions
#-----------------------------------------
region_flags <- c("Northwestern", "Northeastern", "Southern", "Western", "Eastern")
Ken_indices  <- lapply(region_flags[1:3], function(f) which(Ken.Site$flag == f))
Tan_indices  <- lapply(region_flags[4:5], function(f) which(Tan.Site$flag == f))
Score_Data   <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
#-----------------------------------------
#-----------------------------------------
Y_vars <- c("Prevalence")
#-----------------------------------------
#-----------------------------------------
Site.base      <- Site
Ken.Site.base  <- Ken.Site
Tan.Site.base  <- Tan.Site

YEAR      <- unique(Score_Data$Year) %>% sort()
date.time <- data.frame(time.index = 1:length(YEAR),
                        time.scale = seq(0, 1, , length(YEAR)),
                        Year       = YEAR)

Data.base      <- Score_Data %>% left_join(date.time, by = c("Year"))
Data.base$true.Prevalence <- Data.base$Prevalence
cv.Pred.details <- CV.error <- NULL
for(i in 1:50){
  set.seed(i)
  Score_CV_Data <- Data.base
  Site          <- Site.base
  Ken.Site      <- Ken.Site.base
  Tan.Site      <- Tan.Site.base
  setDF(Site)


  Test <- NULL
  for(r in 1:length(as.character(unique(Site$flag)))){
    Village_ID <- Site[(Site$flag %in% unique(Site$flag)[r]) &(Site$Study_Arm %in% c(1, 2, 4)), ]$Village_ID
    Test <- c(Test, sample(Village_ID, 0.25*nrow(Site)*length(Village_ID)/nrow(Site), replace = F))
  }

  setDT(Site)

  Ken.Site$Flag <- ifelse(Ken.Site$Village_ID %in% Test, "test", "train")
  Tan.Site$Flag <- ifelse(Tan.Site$Village_ID %in% Test, "test", "train")


  Score_CV_Data$Flag <- ifelse((Score_CV_Data$Village_ID %in% c(Test)), "test", "train")


  Da.mod  <- Score_CV_Data[Score_CV_Data$Flag == "train", ]
  Da.test <- Score_CV_Data[Score_CV_Data$Flag == "test", ]



  Inla.data <- Score_CV_Data
  Inla.data$Prevalence <- ifelse(Inla.data$Flag == "test", NA, Inla.data$Prevalence)
  pred.ind   <- which(Inla.data$Flag == "test")


  Inla.data$Prevalence <- log(-log(1 - Inla.data$Prevalence))
  setDF(Inla.data)


  nsamp <- 3e2  # which is comparable with an ensemble size of 300
  #JSTVC_xi
  str.1 <- paste0(Y_vars, " ~ -1 + ", paste0(X_vars, collapse = " + "))

  JSTVC_xi  <- inla(as.formula(str.1),
                      family  = "gaussian",
                      data    = Inla.data,
                      verbose = FALSE,
                      control.predictor = list(compute = TRUE),
                      control.compute   = list(config = TRUE))

  post.samp <- inla.posterior.sample(nsamp, JSTVC_xi)

  latent.names <- rownames(post.samp[[1]]$latent)

  # Match predictor indices
  idx <- grep("^Predictor", latent.names)

  # Posterior samples for predictors
  pred.samples <- sapply(
    post.samp,
    function(x) x$latent[idx[pred.ind]]
  )
  pred.orig             <- 1 - exp(-exp(pred.samples)) #exp(5 - exp(-pred.samples))
  Da.test$pred.JSTVC_xi <- apply(pred.orig, c(1), median)
  Da.test$sd.JSTVC_xi   <- apply(pred.orig, c(1), sd)


  #Sub_AR1
  str.2 <- paste0(Y_vars, " ~ -1 + f(Year, model = 'ar1') + f(flag, model = 'iid') + ",
                  paste0(X_vars, collapse = " + "))

  #
  Sub_AR1  <- inla(as.formula(str.2),
                      family = "gaussian",
                      data   = Inla.data,
                      verbose = FALSE,
                      control.predictor = list(compute = TRUE),
                      control.compute   = list(config = TRUE))

  post.samp <- inla.posterior.sample(nsamp, Sub_AR1)

  latent.names <- rownames(post.samp[[1]]$latent)

  # Match predictor indices
  idx <- grep("^Predictor", latent.names)

  # Posterior samples for predictors
  pred.samples <- sapply(
                          post.samp,
                          function(x) x$latent[idx[pred.ind]]
                        )
  pred.orig            <- 1 - exp(-exp(pred.samples))
  Da.test$pred.Sub_AR1 <- apply(pred.orig, c(1), median)
  Da.test$sd.Sub_AR1   <- apply(pred.orig, c(1), sd)



  #arm.AR1
  str.3 <- paste0(Y_vars, " ~ -1 + f(Year, model = 'ar1') + f(Study_Arm, model = 'iid') + ",
                  paste0(X_vars, collapse = " + "))

  #
  Arm_AR1  <- inla(as.formula(str.3),
                      family = "gaussian",
                      data   = Inla.data,
                      verbose = FALSE,
                      control.predictor = list(compute = TRUE),
                      control.compute   = list(config = TRUE))

  post.samp <- inla.posterior.sample(nsamp, Arm_AR1)

  latent.names <- rownames(post.samp[[1]]$latent)

  # Match predictor indices
  idx <- grep("^Predictor", latent.names)

  # Posterior samples for predictors
  pred.samples <- sapply(
    post.samp,
    function(x) x$latent[idx[pred.ind]]
  )
  pred.orig            <- 1 - exp(-exp(pred.samples))
  Da.test$pred.Arm_AR1 <- apply(pred.orig, c(1), median)
  Da.test$sd.Arm_AR1   <- apply(pred.orig, c(1), sd)
  #-----------------------------------------
  # Beta estimation per region
  #-----------------------------------------
  Beta.list    <- vector("list", length(region_flags))
  G.basic.data <- vector("list", length(region_flags))

  for(r in 1:length(region_flags)) {
    if(r <=3){
      Temp <- Ken.Site[Ken_indices[[r]], c("Village_ID","LAT","LON")] %>%
        left_join(Score_CV_Data, by = c("Village_ID"))
    }else{
      Temp <- Tan.Site[Tan_indices[[r - 3]], c("Village_ID","LAT","LON")] %>%
        left_join(Score_CV_Data, by = c("Village_ID"))
    }


    cov.ind <- which(colnames(Temp) %in% c("Village_ID","log.log.mean.Intensity","log.log.var.Intensity"))

    Var.variables <- Temp[, cov.ind]
    Var.variables <- aggregate(as.formula(paste("log(", Y_vars, ")~ Village_ID")), data = Temp, var, na.rm = TRUE) %>%
      left_join(Var.variables, by = "Village_ID") %>%
      unique()

    Var.variables <- Var.variables[, -1]
    colnames(Var.variables) <- c("var.log.log.Prevalence",
                                 "log.log.mean.Intensity",
                                 "log.log.var.Intensity")

    Beta.list[[r]] <- as.vector(coef(lm(var.log.log.Prevalence ~
                                          exp(log.log.mean.Intensity) +
                                          exp(log.log.var.Intensity), data = Var.variables)))

    G.basic.data[[r]] <- regCreateGridm(
      sample.coords   = if(r <= 3) Ken.Site[Ken_indices[[r]], ] else Tan.Site[Tan_indices[[r - 3]], ],
      adjacent.matrix = if(r <= 3) Ken.G else Tan.G,
      loc             = c("LAT","LON"),
      pred.coords     = if(r <= 3) Ken.Site[Ken_indices[[r]], ] else Tan.Site[Tan_indices[[r - 3]], ],
      BAUs.Dist       = if(r <= 3) Kenya.Dist.c[Ken_indices[[r]], Ken_indices[[r]]] else Tanzania.Dist.c[Tan_indices[[r - 3]], Tan_indices[[r - 3]]],
      R.sqrt          = 0,
      site.id         = "Village_ID",
      ch              = 50,
      method          = "Wenland",
      H.Grid_dist     = if(r <= 3) Kenya.Dist.c[Ken_indices[[r]], Ken_indices[[r]]] else Tanzania.Dist.c[Tan_indices[[r - 3]], Tan_indices[[r - 3]]],
      var.covariable  = Var.variables[, -1]
    )
  }
  #-----------------------------------------
  # Prepare fixed effect data per region
  #-----------------------------------------
  Fixed.effect.Data <- lapply(region_flags, function(f) {
    Construct.Fixed.effect.Data(
      data = subset(Score_CV_Data, flag == f),
      date_time = "time.index",
      siteid = "Village_ID",
      include = NULL,
      Y = Y_vars,
      X = X_vars,
      R = NULL,
      Z = NULL,
      G = NULL,
      sX = NULL,
      sR = NULL,
      sZ = NULL,
      sG = c("Intercept"),
      standard = FALSE,
      center = FALSE,
      start.index = 1
    )
  })
  names(Fixed.effect.Data) <- region_flags
  #-----------------------------------------
  # Priors and initial parameters
  #-----------------------------------------
  res.num <- G.basic.data[[1]]$Grid.infor$summary$res
  Py <- length(Fixed.effect.Data)
  Px   <- sapply(Fixed.effect.Data, function(x) if(is.null(x$X_ts)) 0 else nrow(x$X_ts))

  #-----------------------
  # 3. Prior
  #-----------------------
  prior <- lapply(seq_len(Py), function(py){
    list(
      beta  = list(mu.beta  = rep(0, length(X_vars)),
                   sigma.sq = diag(Px[py]) * c(1e4, rep(1e4, Px[py] - 1))),
      obs.sigma.sq = list(a = 2, b = 1),
      st.sRF = list(
        Initial = list(
          rho.v = list(mu = rep(1e-3, res.num), sigma.sq = rep(1e5, res.num)),
          Phi.v = list(a = rep(2, res.num), b = rep(1, res.num)),
          LR.coef = list(mu.delta = Beta.list[[py]][-1], sigma.sq = diag(2)*1e4),
          proc.tau.sq = list(a = rep(2, res.num), b = rep(1, res.num))
        )
      )
    )
  })
  names(prior) <- names(Fixed.effect.Data)
  #-----------------------
  # 4. Para.List
  #-----------------------
  Para.List <- lapply(seq_len(Py), function(py){
    list(
      obs.sigma.sq = list(mu.sigma.sq = 1, a = 2, b = 1),
      st.sRF = list(
        Initial = list(
          rho.v = list(mu.rho.v = rep(1e-1, res.num), var.rho.v = rep(1e-2, res.num)),
          Phi.v = list(mu.Phi.v = 50,
                       L = rep(1e1, res.num),
                       H = rep(1e3, res.num)),
          LR.coef = list(mu.alpha = Beta.list[[py]][-1]),
          proc.tau.sq = list(mu.tau.sq = rep(1e-3, res.num),
                             L = rep(1e-5, res.num),
                             H = rep(1e2, res.num))
        )
      )
    )
  })
  names(Para.List) <- names(Fixed.effect.Data)
  #-----------------------------------------
  # Model fitting
  #-----------------------------------------
  start.time <- proc.time()
  CV_Ranalysis <- JSTVC(Fixed.effect.Data = Fixed.effect.Data,
                        G.basic.data      = G.basic.data,
                        prior             = prior,
                        Para.List         = Para.List,
                        CV                = TRUE,
                        Object            = "Flag",
                        transf.Response   = "loglog",
                        plot              = TRUE,
                        Ne                = 300,
                        tol.real          = 1e-5,
                        itMin             = 1e0,
                        itMax             = 1e3,
                        Obj.Seq           = 1)
  end.time <- proc.time()
  print(end.time - start.time)
  JSTVC.Pred.Res <- JSTVC.Pred.sd <- NULL
  for(py in 1:5){
    Pred.mu <- apply(CV_Ranalysis[[1]]$test.Pred.ens[[py]], c(1, 2), median)
    Pred.sd <- apply(CV_Ranalysis[[1]]$test.Pred.ens[[py]], c(1, 2), sd)

    sRF.1      <- as.data.frame(Pred.mu)
    Col.Name   <- colnames(sRF.1)
    sRF.1$Year <- c(2011, 2012, 2013, 2014, 2015)
    sRF.1 <- sRF.1 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='pred.JSTVC') %>% as.data.frame()

    Pred.sd      <- as.data.frame(Pred.sd)
    Col.Name     <- colnames(Pred.sd)
    Pred.sd$Year <- c(2011, 2012, 2013, 2014, 2015)
    Pred.sd      <- Pred.sd %>% pivot_longer(cols = Col.Name,
                                        names_to ='Village_ID',
                                        values_to ='sd.JSTVC') %>% as.data.frame()


    JSTVC.Pred.Res <- rbind(JSTVC.Pred.Res, sRF.1[(sRF.1$Village_ID %in% c(Test)), ])
    JSTVC.Pred.sd <- rbind(JSTVC.Pred.sd, Pred.sd[(Pred.sd$Village_ID %in% c(Test)), ])
  }


  temp <- JSTVC.Pred.Res %>%
          left_join(JSTVC.Pred.sd, by = c("Village_ID", "Year")) %>%
          left_join(Da.test[, c(1, 2, 45:52)], by = c("Village_ID", "Year"))
  temp$fold <- i
  HDCM::spT_validation()


  library(HDCM)
  JSTVC.spT <- spT_validation(z       = temp$true.Prevalence,
                             zhat     = temp$pred.JSTVC,
                             zhat.Ens = NULL,
                             sigma    = temp$sd.JSTVC,
                             names    = F,
                             CC       = F)

  JSTVC_xi.spT <- spT_validation(z        = temp$true.Prevalence,
                               zhat     = temp$pred.JSTVC_xi,
                               zhat.Ens = NULL,
                               sigma    = temp$sd.JSTVC_xi,
                               names    = F,
                               CC       = F)

  Sub_AR1.spT <- spT_validation(z        = temp$true.Prevalence,
                               zhat     = temp$pred.Sub_AR1,
                               zhat.Ens = NULL,
                               sigma    = temp$sd.Sub_AR1,
                               names    = F,
                               CC       = F)

  Arm_AR1.spT <- spT_validation(z        = temp$true.Prevalence,
                               zhat     = temp$pred.Arm_AR1,
                               zhat.Ens = NULL,
                               sigma    = temp$sd.Arm_AR1,
                               names    = F,
                               CC       = F)
  if(i > 1){
    load(file = Tab.Name)
  }
  CV.error <- rbind(CV.error, data.frame(iter        = i,
                                         STVC.RMSE     = JSTVC.spT[1],
                                         JSTVC_xi.RMSE = JSTVC_xi.spT[1],
                                         Sub_AR1.RMSE  = Sub_AR1.spT[1],
                                         Arm_AR1.RMSE  = Arm_AR1.spT[1],
                                         STVC.CRPS     = JSTVC.spT[4],
                                         JSTVC_xi.CRPS = JSTVC_xi.spT[4],
                                         Sub_AR1.CRPS  = Sub_AR1.spT[4],
                                         Arm_AR1.CRPS  = Arm_AR1.spT[4]

  ))
  cv.Pred.details <- rbind(cv.Pred.details, temp)
  save(CV.error, cv.Pred.details, file = Tab.Name)

  print(CV.error)
}



#-----------------------------------------
#2. JSTVC without IE
#-----------------------------------------
X_vars   <- c("Intercept", paste0("CWT_", 1:4), paste0("SBT_", 1:4))
Tab.Name <- paste0(Tab, "/SCORE_CV_models_RMSE_CRPS_no_ie.RData")
cv.Pred.details <- CV.error <- NULL
for(i in 1:50){
  set.seed(i)
  Score_CV_Data <- Data.base
  Site          <- Site.base
  Ken.Site      <- Ken.Site.base
  Tan.Site      <- Tan.Site.base
  setDF(Site)


  Test <- NULL
  for(r in 1:length(as.character(unique(Site$flag)))){
    Village_ID <- Site[(Site$flag %in% unique(Site$flag)[r]) &(Site$Study_Arm %in% c(1, 2, 4)), ]$Village_ID
    Test <- c(Test, sample(Village_ID, 0.25*nrow(Site)*length(Village_ID)/nrow(Site), replace = F))
  }

  setDT(Site)

  Ken.Site$Flag <- ifelse(Ken.Site$Village_ID %in% Test, "test", "train")
  Tan.Site$Flag <- ifelse(Tan.Site$Village_ID %in% Test, "test", "train")


  Score_CV_Data$Flag <- ifelse((Score_CV_Data$Village_ID %in% c(Test)), "test", "train")
  Da.test <- Score_CV_Data[Score_CV_Data$Flag == "test", ]
  #-----------------------------------------
  # Beta estimation per region
  #-----------------------------------------
  Beta.list    <- vector("list", length(region_flags))
  G.basic.data <- vector("list", length(region_flags))

  for(r in 1:length(region_flags)) {
    if(r <=3){
      Temp <- Ken.Site[Ken_indices[[r]], c("Village_ID","LAT","LON")] %>%
        left_join(Score_CV_Data, by = c("Village_ID"))
    }else{
      Temp <- Tan.Site[Tan_indices[[r - 3]], c("Village_ID","LAT","LON")] %>%
        left_join(Score_CV_Data, by = c("Village_ID"))
    }


    cov.ind <- which(colnames(Temp) %in% c("Village_ID","log.log.mean.Intensity","log.log.var.Intensity"))

    Var.variables <- Temp[, cov.ind]
    Var.variables <- aggregate(as.formula(paste("log(", Y_vars, ")~ Village_ID")), data = Temp, var, na.rm = TRUE) %>%
      left_join(Var.variables, by = "Village_ID") %>%
      unique()

    Var.variables <- Var.variables[, -1]
    colnames(Var.variables) <- c("var.log.log.Prevalence",
                                 "log.log.mean.Intensity",
                                 "log.log.var.Intensity")

    Beta.list[[r]] <- as.vector(coef(lm(var.log.log.Prevalence ~
                                          exp(log.log.mean.Intensity) +
                                          exp(log.log.var.Intensity), data = Var.variables)))

    G.basic.data[[r]] <- regCreateGridm(
      sample.coords   = if(r <= 3) Ken.Site[Ken_indices[[r]], ] else Tan.Site[Tan_indices[[r - 3]], ],
      adjacent.matrix = if(r <= 3) Ken.G else Tan.G,
      loc             = c("LAT","LON"),
      pred.coords     = if(r <= 3) Ken.Site[Ken_indices[[r]], ] else Tan.Site[Tan_indices[[r - 3]], ],
      BAUs.Dist       = if(r <= 3) Kenya.Dist.c[Ken_indices[[r]], Ken_indices[[r]]] else Tanzania.Dist.c[Tan_indices[[r - 3]], Tan_indices[[r - 3]]],
      R.sqrt          = 0,
      site.id         = "Village_ID",
      ch              = 50,
      method          = "Wenland",
      H.Grid_dist     = if(r <= 3) Kenya.Dist.c[Ken_indices[[r]], Ken_indices[[r]]] else Tanzania.Dist.c[Tan_indices[[r - 3]], Tan_indices[[r - 3]]],
      var.covariable  = Var.variables[, -1]
    )
  }
  #-----------------------------------------
  # Prepare fixed effect data per region
  #-----------------------------------------
  Fixed.effect.Data <- lapply(region_flags, function(f) {
    Construct.Fixed.effect.Data(
      data = subset(Score_CV_Data, flag == f),
      date_time = "time.index",
      siteid = "Village_ID",
      include = NULL,
      Y = Y_vars,
      X = X_vars,
      R = NULL,
      Z = NULL,
      G = NULL,
      sX = NULL,
      sR = NULL,
      sZ = NULL,
      sG = c("Intercept"),
      standard = FALSE,
      center = FALSE,
      start.index = 1
    )
  })
  names(Fixed.effect.Data) <- region_flags
  #-----------------------------------------
  # Priors and initial parameters
  #-----------------------------------------
  res.num <- G.basic.data[[1]]$Grid.infor$summary$res
  Py <- length(Fixed.effect.Data)
  Px   <- sapply(Fixed.effect.Data, function(x) if(is.null(x$X_ts)) 0 else nrow(x$X_ts))

  #-----------------------
  # 3. Prior
  #-----------------------
  prior <- lapply(seq_len(Py), function(py){
    list(
      beta  = list(mu.beta  = rep(0, length(X_vars)),
                   sigma.sq = diag(Px[py]) * c(1e4, rep(1e4, Px[py] - 1))),
      obs.sigma.sq = list(a = 2, b = 1),
      st.sRF = list(
        Initial = list(
          rho.v = list(mu = rep(1e-3, res.num), sigma.sq = rep(1e5, res.num)),
          Phi.v = list(a = rep(2, res.num), b = rep(1, res.num)),
          LR.coef = list(mu.delta = Beta.list[[py]][-1], sigma.sq = diag(2)*1e4),
          proc.tau.sq = list(a = rep(2, res.num), b = rep(1, res.num))
        )
      )
    )
  })
  names(prior) <- names(Fixed.effect.Data)
  #-----------------------
  # 4. Para.List
  #-----------------------
  Para.List <- lapply(seq_len(Py), function(py){
    list(
      obs.sigma.sq = list(mu.sigma.sq = 1, a = 2, b = 1),
      st.sRF = list(
        Initial = list(
          rho.v = list(mu.rho.v = rep(1e-1, res.num), var.rho.v = rep(1e-2, res.num)),
          Phi.v = list(mu.Phi.v = 50,
                       L = rep(1e1, res.num),
                       H = rep(1e3, res.num)),
          LR.coef = list(mu.alpha = Beta.list[[py]][-1]),
          proc.tau.sq = list(mu.tau.sq = rep(1e-3, res.num),
                             L = rep(1e-5, res.num),
                             H = rep(1e2, res.num))
        )
      )
    )
  })
  names(Para.List) <- names(Fixed.effect.Data)
  #-----------------------------------------
  # Model fitting
  #-----------------------------------------
  start.time <- proc.time()
  CV_Ranalysis <- JSTVC(Fixed.effect.Data = Fixed.effect.Data,
                        G.basic.data      = G.basic.data,
                        prior             = prior,
                        Para.List         = Para.List,
                        CV                = TRUE,
                        Object            = "Flag",
                        transf.Response   = "loglog",
                        plot              = TRUE,
                        Ne                = 300,
                        tol.real          = 1e-5,
                        itMin             = 1e0,
                        itMax             = 1e3,
                        Obj.Seq           = 1)
  end.time <- proc.time()
  print(end.time - start.time)
  JSTVC.Pred.Res <- JSTVC.Pred.sd <- NULL
  for(py in 1:5){
    Pred.mu <- apply(CV_Ranalysis[[1]]$test.Pred.ens[[py]], c(1, 2), median)
    Pred.sd <- apply(CV_Ranalysis[[1]]$test.Pred.ens[[py]], c(1, 2), sd)

    sRF.1      <- as.data.frame(Pred.mu)
    Col.Name   <- colnames(sRF.1)
    sRF.1$Year <- c(2011, 2012, 2013, 2014, 2015)
    sRF.1 <- sRF.1 %>% pivot_longer(cols = Col.Name,
                                    names_to ='Village_ID',
                                    values_to ='pred.JSTVC') %>% as.data.frame()

    Pred.sd      <- as.data.frame(Pred.sd)
    Col.Name     <- colnames(Pred.sd)
    Pred.sd$Year <- c(2011, 2012, 2013, 2014, 2015)
    Pred.sd      <- Pred.sd %>% pivot_longer(cols = Col.Name,
                                             names_to ='Village_ID',
                                             values_to ='sd.JSTVC') %>% as.data.frame()


    JSTVC.Pred.Res <- rbind(JSTVC.Pred.Res, sRF.1[(sRF.1$Village_ID %in% c(Test)), ])
    JSTVC.Pred.sd <- rbind(JSTVC.Pred.sd, Pred.sd[(Pred.sd$Village_ID %in% c(Test)), ])
  }


  temp <- JSTVC.Pred.Res %>%
    left_join(JSTVC.Pred.sd, by = c("Village_ID", "Year")) %>%
    left_join(Da.test[, c(1, 2, 45:46)], by = c("Village_ID", "Year"))
  temp$fold <- i

  JSTVC_ie.spT <- spT_validation(z       = temp$true.Prevalence,
                              zhat     = temp$pred.JSTVC,
                              zhat.Ens = NULL,
                              sigma    = temp$sd.JSTVC,
                              names    = F,
                              CC       = F)
  if(i > 1){
    load(file = Tab.Name)
  }
  CV.error <- rbind(CV.error, data.frame(iter         = i,
                                         STVC_ie.RMSE = JSTVC_ie.spT[1],
                                         STVC_ie.CRPS = JSTVC_ie.spT[4]))
  cv.Pred.details <- rbind(cv.Pred.details, temp)
  save(CV.error, cv.Pred.details, file = Tab.Name)
  print(CV.error)
}
