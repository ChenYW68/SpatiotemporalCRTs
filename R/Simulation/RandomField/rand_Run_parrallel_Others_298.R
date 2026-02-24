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
  "spam",
  "INLA"
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
start  <- c(1, 300)
n.cores <- 15
cl      <- makeCluster(n.cores)
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
                     "Tanzania_Score_Data"
))

clusterEvalQ(cl, {
  source(normalizePath("./JSTVC/R/util.R"))
  source(normalizePath("./JSTVC/R/sim_Generate_Data.R"))
})

Tab <- paste0("./result/Simulation_300/random_competing_n_", 298)
if (!dir.exists(Tab)) {
  dir.create(Tab, recursive = TRUE)
}
#-----------------------------------------
# Parrallel
#-----------------------------------------

clusterExport(cl, c("Tab"))
clusterExport(cl, ls())
#-----------------------------------------
results_list <- parLapply(cl, start[1]:start[2], function(iter) {
  start.time <- proc.time()
  set.seed(iter)

  # -------------------- Main --------------------
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

  simData <- Simu_stData(para = sim.para,
                         loc  = Simu_data$loc,
                         W_ts = Simu_data$W_ts,
                         X    = sim.para$X)
  sim_Data         <- data.collect(simData, sim.para)
  Train.village.ID <- unique(sim_Data$Village_ID[sim_Data$Simu == "Train"])
  sim_Data$n.size  <- length(Train.village.ID)


  dummies           <- model.matrix(~ factor(flag) - 1, data = sim_Data)
  colnames(dummies) <- paste0("X", 1:ncol(dummies))
  sim_Data          <- cbind(sim_Data, dummies)

  # train/test flags
  Ken.Site$Flag <- ifelse(Ken.Site$Village_ID %in% Train.village.ID, "train", "test")
  Tan.Site$Flag <- ifelse(Tan.Site$Village_ID %in% Train.village.ID, "train", "test")
  sim_Data$Flag <- ifelse(sim_Data$Village_ID %in% Train.village.ID, "train", "test")


  Da.mod <- sim_Data %>% left_join(Score_Data[, c(1, 2, 40:41)], by = c("Village_ID", "Year"))


  X_vars     <- c("Intercept", colnames(sim.para$X)[-1])
  str.1 <- paste0("Y_ts ~ -1 + ", paste0(X_vars, collapse = " + "))

  str.2 <- paste0("Y_ts ~ -1 + f(Year, model = 'ar1') + f(flag, model = 'iid') + ",
                  paste0(X_vars, collapse = " + "))

  str.3 <- paste0("Y_ts ~ -1 + f(Year, model = 'ar1') + f(Study_Arm, model = 'iid') + ",
                  paste0(X_vars, collapse = " + "))


  JSTVC.nonEffect  <- inla(as.formula(str.1),
                           family  = "gaussian",
                           data    = Da.mod,
                           verbose = FALSE,
                           control.predictor = list(compute = TRUE))
  JSTVC.subRegion  <- inla(as.formula(str.2),
                           family  = "gaussian",
                           data    = Da.mod,
                           verbose = FALSE,
                           control.predictor = list(compute = TRUE))
  JSTVC.subArm  <- inla(as.formula(str.3),
                        family  = "gaussian",
                        data    = Da.mod,
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
       file = paste0(Tab, "/sim_", iter, ".RData"))
  return(1)
})
stopCluster(cl)
