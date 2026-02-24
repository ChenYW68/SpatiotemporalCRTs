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

start  <- c(1, 300)   # 
px     <- 11
Y_vars <- c("Y_ts")
X_vars <- c("Intercept", paste0("CWT_", 1:4), paste0("SBT_", 1:4),
            "spillover.CWT", "spillover.SBT")

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
  Tab <- paste0("./Result/Simulation/VB/smoothed_other_n_", n,
                "_Ch_", Ch[1],
                "_Cs_", Cs[1],
                "_Ne_", Ne)

  if (!dir.exists(Tab)) {
    dir.create(Tab, recursive = TRUE)
  }
  #-----------------------------------------
  #
  #-----------------------------------------

  clusterExport(cl, c(
    "Tab", "Cs", "Ne"))
  clusterExport(cl, ls())
  #-----------------------------------------
  # 
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



    # Da.mod <- sim_Data#[sim_Data$Flag %in% "train",]


    X_vars     <- c("Intercept", colnames(sim.para$X)[-1])
    str.1 <- paste0("Y_ts ~ -1 + ", paste0(X_vars, collapse = " + "))

    str.2 <- paste0("Y_ts ~ -1 + f(Year, model = 'ar1') + f(flag, model = 'iid') + ",
                    paste0(X_vars, collapse = " + "))

    str.3 <- paste0("Y_ts ~ -1 + f(Year, model = 'ar1') + f(Study_Arm, model = 'iid') + ",
                    paste0(X_vars, collapse = " + "))


    JSTVC.nonEffect  <- inla(as.formula(str.1),
                       family  = "gaussian",
                       data    = sim_Data,
                       verbose = FALSE,
                       control.predictor = list(compute = TRUE))
    JSTVC.subRegion  <- inla(as.formula(str.2),
                             family  = "gaussian",
                             data    = sim_Data,
                             verbose = FALSE,
                             control.predictor = list(compute = TRUE))
    JSTVC.subArm  <- inla(as.formula(str.3),
                          family  = "gaussian",
                          data    = sim_Data,
                          verbose = FALSE,
                          control.predictor = list(compute = TRUE))

    JSTVC.nonEffect.beta <- round(summary(JSTVC.nonEffect)[[3]][c(1:length(X_vars)), 1], 3)
    JSTVC.nonEffect.sd <- round(summary(JSTVC.nonEffect)[[3]][c(1:length(X_vars)), 2], 3)
    names(JSTVC.nonEffect.beta) <-  names(JSTVC.nonEffect.sd) <- paste0("nonEffect.", X_vars)

    JSTVC.subRegion.beta <- round(summary(JSTVC.subRegion)[[3]][c(1:length(X_vars)), 1], 3)
    JSTVC.subRegion.sd <- round(summary(JSTVC.subRegion)[[3]][c(1:length(X_vars)), 2], 3)
    names(JSTVC.subRegion.beta) <- names(JSTVC.subRegion.sd) <- paste0("subRegion.", X_vars)

    JSTVC.subArm.beta <- round(summary(JSTVC.subArm)[[3]][c(1:length(X_vars)), 1], 3)
    JSTVC.subArm.sd <- round(summary(JSTVC.subArm)[[3]][c(1:length(X_vars)), 2], 3)
    names(JSTVC.subArm.beta) <- names(JSTVC.subArm.sd) <- paste0("subArm.", X_vars)




    save(JSTVC.nonEffect.beta, JSTVC.nonEffect.sd,
    JSTVC.subRegion.beta, JSTVC.subRegion.sd,
    JSTVC.subArm.beta, JSTVC.subArm.sd,
    file = paste0(Tab, "/sim_", iter, "_other.RData"))
    return(1)
  })
}
stopCluster(cl)
