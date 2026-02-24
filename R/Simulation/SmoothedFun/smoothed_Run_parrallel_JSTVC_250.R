#-----------------------------------------
# Clear workspace and load environment
#-----------------------------------------
rm(list = ls())
library(parallel)
library(data.table)
library(dplyr)
source("./LoadPackages/RDependPackages.R")
# Packages to load
pkgs <- c(
  "data.table",
  "ggplot2",
  "plyr",
  "parallel",
  "sqldf",
  "numDeriv",
  "lubridate",
  "dplyr",
  "Hmisc",
  "MASS",
  "tidyr",
  "RColorBrewer",
  "progress",
  "fields",
  "ranger",
  "MBA",
  "Rcpp",
  "writexl",
  "readxl",
  "ggmap",
  "verification",
  "mapproj",
  "sp",
  "mvnfast",
  "INLA",
  "spam",
  "scoringutils",
  "HDCM",
  "rARPACK"
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
n   <- 250
Ch  <- rep(50, 5)

start  <- c(1, 300) 
px     <- 11
Y_vars <- c("Y_ts")
X_vars <- c("Intercept"
            , paste0("CWT_", 1:4)
            , paste0("SBT_", 1:4)
            , "IEt.CWT"
            , "IEt.SBT"
)

grid <- rbind(expand.grid(
  # ch = seq(20, 100, by = 30),
  # cs = 50,#seq(20, 50, by = 5), #50, 300
  # cs = seq(50, 300, by = 50),
  ne = seq(50, 500, by = 10)
)
# , data.frame(cs = c(150, 200, 250, 300),
#              ne = c(500, 500, 500, 500))
)



n.cores <- 15
# n.cores <- parallel::detectCores() - 1

cl      <- makeCluster(n.cores)
# Assuming you have a cluster created, e.g.,
# cl <- makeCluster(detectCores())
clusterExport(cl, "pkgs")
# Load all packages on all cluster nodes
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
                     "n",
                     "Ch",
                     # "Cs",
                     # "Ne",
                     "px",
                     "Y_vars",
                     "X_vars",
                     "Kenya_Score_Data",
                     "Tanzania_Score_Data"
))

clusterEvalQ(cl, {
  source(normalizePath("./nSTJVC/R/regCreateGridm2.R"))
  source(normalizePath("./nSTJVC/R/Partitioning.Datasets.R"))
  source(normalizePath("./nSTJVC/R/Construct.Fixed.effect.Data.R"))
  source(normalizePath("./nSTJVC/R/util_08_10.R"))

  source(normalizePath("./nSTJVC/R/VB.Laplace_08_10.R"))
  # source(normalizePath("./nSTJVC/MCMC/MCMC.R"))


  source(normalizePath("./nSTJVC/R/VB.LA.spAugEnKS.R"))
  source(normalizePath("./nSTJVC/R/spAugEnKS.R"))
  source(normalizePath("./nSTJVC/R/JSTVCs.R"))
  source(normalizePath("./simulation/sim_Generate_Data.R"))
})

# stopCluster(cl)#nrow(grid)
for(cv in 1:1){
  Cs  <- rep(1e15, 5) #grid$cs[cv]
  Ne  <- 300#grid$ne[cv]
  Tab <- paste0("./Result/Simulation/VB/smoothed_vb_n_", n,
                "_Ch_", Ch[1],
                "_Cs_", Cs[1],
                "_Ne_", Ne)

  if (!dir.exists(Tab)) {
    dir.create(Tab, recursive = TRUE)
  }

  clusterExport(cl, c(
    "Tab", "Cs", "Ne"))
  clusterExport(cl, ls())
  #-----------------------------------------

  #-----------------------------------------
  results_list <- parLapply(cl, start[1]:start[2], function(iter) {
    # source(normalizePath("./LoadPackages/RDependPackages.R"))
    # source(normalizePath("./nSTJVC/R/regCreateGridm2.R"))
    # source(normalizePath("./nSTJVC/R/Partitioning.Datasets.R"))
    # source(normalizePath("./nSTJVC/R/Construct.Fixed.effect.Data.R"))
    # source(normalizePath("./nSTJVC/R/util_08_10.R"))
    #
    # source(normalizePath("./nSTJVC/R/VB.Laplace_08_10.R"))
    # # source(normalizePath("./nSTJVC/MCMC/MCMC.R"))
    #
    #
    # source(normalizePath("./nSTJVC/R/VB.LA.spAugEnKS.R"))
    # source(normalizePath("./nSTJVC/R/spAugEnKS.R"))
    # source(normalizePath("./nSTJVC/R/JSTVCs.R"))
    # source(normalizePath("./simulation/sim_Generate_Data.R"))
    # py <- 1

    start.time <- proc.time()
    set.seed(iter)

    # --------------------  --------------------
    para             <- list(Nt = 5, nugget = 0)
    Simu_data        <- Simu_stData.surface(para)


    # simData.DataBase <- Simu_data$sim.Wts
    Score_Data       <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
    ind.x            <- unlist(lapply(region_flags, function(f) which(Score_Data$flag == f)))

    sim.para <- list(
      n      = 250,
      Nt     = 5,
      px     = px,
      alpha  = c(5, rep(-1, px - 1)),
      nugget = 1e-2,
      X      = Score_Data[ind.x, c(2, 9:16, 19:20)]
    )

    simData <- Simu_stData.surface(para = sim.para,
                                   loc  = Simu_data$loc,
                                   W_ts = Simu_data$W_ts - mean(Simu_data$W_ts),
                                   X    = sim.para$X)
    sim_Data         <- data.collect(simData, sim.para)
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
    sim_Data <- sim_Data %>%
      left_join(Score_Data[, c(1,2, 40:41)], by = c("Village_ID", "Year"))





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

      G.basic.data[[r]] <- regCreateGridm2(
        sample.coords   = if(r <= 3) Ken.Site[Ken_indices[[r]], ] else Tan.Site[Tan_indices[[r - 3]], ],
        adjacent.matrix = if(r <= 3) Ken.G else Tan.G,
        loc             = c("LAT","LON"),
        pred.coords     = if(r <= 3) Ken.Site[Ken_indices[[r]], ] else Tan.Site[Tan_indices[[r - 3]], ],
        BAUs.Dist       = if(r <= 3) Kenya.Dist.c[Ken_indices[[r]], Ken_indices[[r]]] else Tanzania.Dist.c[Tan_indices[[r - 3]], Tan_indices[[r - 3]]],
        R.sqrt          = 0,
        site.id         = "Village_ID",
        ch              = Ch[r],
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
        Y = Y_vars,
        X = X_vars,
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
    #-----------------------
    # 1. 
    #-----------------------
    res.num <- G.basic.data[[1]]$Grid.infor$summary$res
    Py <- length(Fixed.effect.Data)
    # Phi.v <- rep(5e1, 5)

    #-----------------------
    # 2. 
    #-----------------------
    Px   <- sapply(Fixed.effect.Data, function(x) if(is.null(x$X_ts)) 0 else nrow(x$X_ts))
    Pz   <- sapply(Fixed.effect.Data, function(x) if(is.null(x$Z_ts)) 0 else nrow(x$Z_ts))
    Pg   <- sapply(Fixed.effect.Data, function(x) if(is.null(x$G_ts)) 0 else nrow(x$G_ts))
    sPx  <- sapply(Fixed.effect.Data, function(x) if(is.null(x$sX_ts)) 0 else nrow(x$sX_ts))
    sPz  <- sapply(Fixed.effect.Data, function(x) if(is.null(x$sZ_ts)) 0 else nrow(x$sZ_ts))
    sPg  <- sapply(Fixed.effect.Data, function(x) if(is.null(x$sG_ts)) 0 else nrow(x$sG_ts))

    #-----------------------
    # 3. Prior 
    #-----------------------
    prior <- lapply(seq_len(Py), function(py){
      list(
        beta  = list(mu.beta = rep(0, length(X_vars)),
                     sigma.sq = diag(Px[py]) * c(1e4, rep(1e4, Px[py] - 1)),
                     a = rep(1, Px[py]), b = rep(1, Px[py])),
        sbeta = list(mu.beta = rep(50, sPx[py]),
                     sigma.sq = diag(1, sPx[py]),
                     a = rep(1, sPx[py]), b = rep(1, sPx[py])),
        obs.sigma.sq = list(a = 2, b = 1),
        rho.alpha    = list(mu = matrix(1, nrow = Pz[py], ncol = 1), var = matrix(1e-1, nrow = Pz[py], ncol = 1)),
        alpha.tau.sq = list(mu = matrix(1e-1, nrow = Pz[py], ncol = 1), var = matrix(1e-1, nrow = Pz[py], ncol = 1)),
        ini.alpha    = list(mu = matrix(1e-1, nrow = Pz[py], ncol = 1), var = matrix(1e-1, nrow = Pz[py], ncol = 1)),
        ini.alpha.tau.sq = list(mu = matrix(1e-1, nrow = Pz[py], ncol = 1), var = matrix(1e-1, nrow = Pz[py], ncol = 1)),
        self.rho.alpha = list(mu = matrix(1, nrow = sPz[py], ncol = 1), var = matrix(1e-1, nrow = sPz[py], ncol = 1)),
        self.alpha.tau.sq = list(mu = matrix(1e-1, nrow = sPz[py], ncol = 1), var = matrix(1e-1, nrow = sPz[py], ncol = 1)),
        self.ini.alpha   = list(mu = matrix(1e-1, nrow = sPz[py], ncol = 1), var = matrix(1e-1, nrow = sPz[py], ncol = 1)),
        self.ini.alpha.tau.sq = list(mu = matrix(1e-1, nrow = sPz[py], ncol = 1), var = matrix(1e-1, nrow = sPz[py], ncol = 1)),
        st.sRF = list(
          Initial = list(
            rho.v = list(mu = rep(0, res.num), sigma.sq = rep(1e5, res.num)),
            Phi.v = list(a = rep(2, res.num), b = rep(1, res.num)),
            LR.coef = list(mu.delta = Beta.list[[py]][-1], sigma.sq = diag(2)*1e4),
            proc.tau.sq = list(a = rep(2, res.num), b = rep(1, res.num)),
            proc.ini.tau.sq = list(a = rep(2, res.num), b = rep(1, res.num))
          )
        )
      )
    })
    names(prior) <- names(Fixed.effect.Data)

    #-----------------------
    # 4. Para.List
    #-----------------------
    Para.List <- lapply(seq_len(Py), function(py){
      ini.sigma.sq <- 1e0
      list(
        beta  = list(a = rep(2, Px[py]), b = rep(1, Px[py]), p1 = rep(0.5, Px[py])),
        obs.sigma.sq = list(mu.sigma.sq = ini.sigma.sq, a = 2, b = 1),
        rho.alpha = list(mu.rho.alpha = matrix(1, nrow = Pz[py], ncol = 1)),
        alpha.tau.sq = list(mu.alpha.tau.sq = matrix(1, nrow = Pz[py], ncol = 1)),
        ini.alpha = list(mu.ini.alpha = matrix(0, nrow = Pz[py], ncol = 1)),
        ini.alpha.tau.sq = list(mu.ini.alpha.tau.sq = matrix(1, nrow = Pz[py], ncol = 1)),
        self.rho.alpha = list(mu.rho.alpha = matrix(1, nrow = sPz[py], ncol = 1)),
        self.alpha.tau.sq = list(mu.alpha.tau.sq = matrix(1, nrow = sPz[py], ncol = 1)),
        self.ini.alpha = list(mu.ini.alpha = matrix(0, nrow = sPz[py], ncol = 1)),
        self.ini.alpha.tau.sq = list(mu.ini.alpha.tau.sq = matrix(1, nrow = sPz[py], ncol = 1)),
        st.sRF = list(
          Initial = list(
            rho.v = list(mu.rho.v = rep(1e-1, res.num), var.rho.v = rep(1e-2, res.num)),
            Phi.v = list(mu.Phi.v = Ch[py], #a = 2, b = 1,
                         L = rep(1e0, res.num),
                         H = rep(1e3, res.num)),
            LR.coef = list(mu.alpha = Beta.list[[py]][-1]),
            proc.tau.sq = list(mu.tau.sq = rep(1e-3, res.num),
                               L = rep(1e-10, res.num),
                               H = rep(1e-2, res.num)),
            proc.ini.tau.sq = list(mu.tau.sq = rep(1e-3, res.num))
          )
        )
      )
    })
    names(Para.List) <- names(Fixed.effect.Data)


    #-----------------------------------------
    # Model fitting
    #-----------------------------------------
    # tab <- paste0("JSTVC_", sum(G.basic.data[[1]]$Grid.infor$summary$Knots.count))

    CV_Ranalysis <- JSTVC(Tab               = "JSTVC",
                          Fixed.effect.Data = Fixed.effect.Data,
                          G.basic.data      = G.basic.data,
                          center.y          = F,
                          scale.y           = F,
                          prior             = prior,
                          Para.List         = Para.List,
                          CV                = TRUE,
                          verbose.VB        = TRUE,
                          verbose           = TRUE,
                          Object            = "Flag",
                          transf.Response   = "normal",#
                          plot              = FALSE,
                          Database          = NULL,
                          save.Predict      = T,
                          Ne                = Ne,
                          n.cores           = 1,
                          cs                = Cs,
                          ct                = 0,
                          tol.real          = 1e-5,
                          itMin             = 1E0,
                          itMax             = 1E3,
                          # nu                = c(1e-5, 1e-5, 1e-5),
                          # var.select        = FALSE,
                          positive          = TRUE,
                          # threshold         = c(0.5, 0.5),
                          Obj.Seq           = 1)

    # #-----------------------------------------
    # temp <- setDF(CV_Ranalysis$Process.monitoring$VB.spEnKS.Iter)
    # # plot(temp[, 2])
    #
    # beta <- CV_Ranalysis$update.Para.List$Northwestern$beta$mu.beta[, 1]
    # sigma.sq <- CV_Ranalysis$update.Para.List$Northwestern$beta$cov.prob$SD^2

    if(is.null(CV_Ranalysis$update.Para.List)){
      beta   <- CV_Ranalysis[[1]]$update.Para.List$Northwestern$beta$mu.beta[, 1]
      sigma.sq <- CV_Ranalysis[[1]]$update.Para.List$Northwestern$beta$cov.prob$SD^2
    }else{
      beta     <- CV_Ranalysis$update.Para.List$Northwestern$beta$mu.beta[, 1]
      sigma.sq <- CV_Ranalysis$update.Para.List$Northwestern$beta$cov.prob$SD^2
    }
    #
    #
    # # Example usage:
    # plot_normal(mean = mu[1], sd = Sd[1])
    #
    #
    # #-----------------------------------------
    # end.time <- proc.time()
    # run_time <- (end.time - start.time)[3]
    # beta     <- CV_Ranalysis$update.Para.List[[1]]$beta$mu.beta
    # list(Result = CV_Ranalysis, iter = iter, run_time = run_time)
    save(beta, sigma.sq, file = paste0(Tab, "/sim_", iter, ".RData"))

    return(1)

  })
}
stopCluster(cl)

