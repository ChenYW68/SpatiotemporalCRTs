#-----------------------------------------
# Clear workspace and load environment
#-----------------------------------------
rm(list = ls())
source("./LoadPackages/RDependPackages.R")
pkgs <- c(
         "data.table",
          "parallel",
          "lubridate",
          "dplyr",
          "Hmisc",
          "tidyr",
          "fields",
          "ranger",
          "mvnfast",
          "spam"
         )
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
#-----------------------------------------
# Simulation parameters
#-----------------------------------------
start   <- c(1, 300)
cl      <- makeCluster(15)
clusterExport(cl, "pkgs")
clusterEvalQ(cl, {
  lapply(pkgs, require, character.only = TRUE)
})
clusterExport(cl, c( "Ken_indices",
                     "Tan_indices",
                     "region_flags",
                     "Ken.Site",
                     "Tan.Site",
                     "Ken.G",
                     "Tan.G",
                     "region_flags",
                     "Kenya_Score_Data",
                     "Tanzania_Score_Data"))

clusterEvalQ(cl, {
  source(normalizePath("./JSTVC/R/regCreateGridm.R"))
  source(normalizePath("./JSTVC/R/Partitioning.Datasets.R"))
  source(normalizePath("./JSTVC/R/Construct.Fixed.effect.Data.R"))
  source(normalizePath("./JSTVC/R/util.R"))
  source(normalizePath("./JSTVC/R/VB.R"))
  source(normalizePath("./JSTVC/R/VB_EnKF.R"))
  source(normalizePath("./JSTVC/R/EnKF.R"))
  source(normalizePath("./JSTVC/R/JSTVCs.R"))
  source(normalizePath("./JSTVC/R/sim_Generate_Data.R"))
})

for(cv in 1:1){
  Tab <- paste0("./result/Simulation_300/random_JSTVC_n_", 298, "_Ne_", 300)
  if (!dir.exists(Tab)) {
    dir.create(Tab, recursive = TRUE)
  }
  clusterExport(cl, c("Tab"))
  clusterExport(cl, ls())
  results_list <- parLapply(cl, start[1]:start[2], function(iter) {
    start.time <- proc.time()

    set.seed(iter)

    # ----------------------------------------
    para             <- list(Nt = 5)
    Simu_data        <- Simu_stData(para)

    Score_Data       <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
    ind.x            <- unlist(lapply(region_flags, function(f) which(Score_Data$flag == f)))

    sim.para <- list(
      n      = 298,
      Nt     = 5,
      px     = 11,
      alpha  = c(5, rep(-1, 10)),
      nugget = 1e-2,
      X      = Score_Data[ind.x, c(2, 9:16, 19:20)]
    )

    simData   <- Simu_stData(para = sim.para,
                             loc  = Simu_data$loc,
                             W_ts = Simu_data$W_ts,
                             X    = sim.para$X)
    sim_Data  <- data.collect(simData, sim.para)

    # Wether use miss-specifications for decay function
    # if(misspecified_decay){
    #   sim_Data$IEt.CWT <- Score_Data[ind.x, ]$IEt.CWT.exp
    #   sim_Data$IEt.SBT <- Score_Data[ind.x, ]$IEt.SBT.exp
    # }

    Train.village.ID <- unique(sim_Data$Village_ID[sim_Data$Simu == "Train"])
    sim_Data$n.size  <- length(Train.village.ID)

    # dummy variables
    dummies           <- model.matrix(~ factor(flag) - 1, data = sim_Data)
    colnames(dummies) <- paste0("X", 1:ncol(dummies))
    sim_Data          <- cbind(sim_Data, dummies)

    # train/test flags
    Ken.Site$Flag <- ifelse(Ken.Site$Village_ID %in% Train.village.ID, "train", "test")
    Tan.Site$Flag <- ifelse(Tan.Site$Village_ID %in% Train.village.ID, "train", "test")
    sim_Data$Flag <- ifelse(sim_Data$Village_ID %in% Train.village.ID, "train", "test")

    # merge covariates
    sim_Data <- sim_Data %>% left_join(Score_Data[, c(1,2, 40:41)], by = c("Village_ID", "Year"))

    #-----------------------------------------
    # Beta estimation per region
    #-----------------------------------------
    Beta.list    <- vector("list", length(region_flags))
    G.basic.data <- vector("list", length(region_flags))

    for(r in 1:length(region_flags)) {
      if(r <=3){
        Temp <- Ken.Site[Ken_indices[[r]], c("Village_ID","LAT","LON")] %>%
          left_join(sim_Data, by = c("Village_ID"))
      }else{
        Temp <- Tan.Site[Tan_indices[[r - 3]], c("Village_ID","LAT","LON")] %>%
          left_join(sim_Data, by = c("Village_ID"))
      }

      cov.ind <- which(colnames(Temp) %in% c("Village_ID","log.log.mean.Intensity","log.log.var.Intensity"))

      Var.variables <- Temp[, cov.ind]
      Var.variables <- aggregate(Y_ts ~ Village_ID, data = Temp, var, na.rm = TRUE) %>%
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
        data = subset(sim_Data, flag == f),
        date_time = "time.index",
        siteid = "Village_ID",
        include = NULL,
        Y = c("Y_ts"),
        X = c("Intercept"
              , paste0("CWT_", 1:4)
              , paste0("SBT_", 1:4)
              , "IEt.CWT"
              , "IEt.SBT"),
        R = NULL,
        Z = NULL,
        G = NULL,
        sX = NULL,
        sR = NULL,
        sZ = NULL,
        sG = c("Intercept"),
        initial.miss = "TRUE.Y_ts",
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
        beta  = list(mu.beta = rep(0, Px[py]),
                     sigma.sq = diag(Px[py]) * c(1e4, rep(1e4, Px[py] - 1))),
        obs.sigma.sq = list(a = 2, b = 1),
        st.sRF = list(
          Initial = list(
            rho.v = list(mu = rep(0, res.num), sigma.sq = rep(1e5, res.num)),
            Phi.v = list(a = rep(5e1, res.num), b = rep(1, res.num)),
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
        obs.sigma.sq = list(mu.sigma.sq = 1, a = 1, b = 1),
        st.sRF = list(
          Initial = list(
            rho.v = list(mu.rho.v = rep(1e-1, res.num), var.rho.v = rep(1e-2, res.num)),
            Phi.v = list(mu.Phi.v = 50,
                         L = rep(1e0, res.num),
                         H = rep(1e3, res.num)),
            LR.coef = list(mu.alpha = Beta.list[[py]][-1]),
            proc.tau.sq = list(mu.tau.sq = rep(1e-3, res.num),
                               L = rep(1e-5, res.num),
                               H = rep(1e-2, res.num))
          )
        )
      )
    })
    names(Para.List) <- names(Fixed.effect.Data)
    #-----------------------------------------
    # Model fitting
    #-----------------------------------------
    CV_Ranalysis <- JSTVC(Fixed.effect.Data = Fixed.effect.Data,
                          G.basic.data      = G.basic.data,
                          prior             = prior,
                          Para.List         = Para.List,
                          CV                = FALSE,
                          Object            = "Flag",
                          transf.Response   = "normal",
                          plot              = TRUE,
                          Ne                = 300,
                          tol.real          = 1e-5,
                          itMin             = 1e2,
                          itMax             = 1e3,
                          Obj.Seq           = 1)

    # #-----------------------------------------
    if(is.null(CV_Ranalysis$update.Para.List)){
      beta   <- CV_Ranalysis[[1]]$update.Para.List$Northwestern$beta$mu.beta[, 1]
      sigma.sq <- CV_Ranalysis[[1]]$update.Para.List$Northwestern$beta$cov.prob$SD^2
    }else{
      beta     <- CV_Ranalysis$update.Para.List$Northwestern$beta$mu.beta[, 1]
      sigma.sq <- CV_Ranalysis$update.Para.List$Northwestern$beta$cov.prob$SD^2
    }
    save(beta, sigma.sq, file = paste0(Tab, "/sim_", iter, ".RData"))
    return(1)
  })
}
stopCluster(cl)
