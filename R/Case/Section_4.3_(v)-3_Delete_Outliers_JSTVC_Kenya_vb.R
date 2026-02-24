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
X_vars <- c("Intercept", paste0("CWT_", 1:4), paste0("SBT_", 1:4), "IEt.CWT", "IEt.SBT")
Tab.Name <- paste0(Tab, "/Delete_Outliers_JSTVC_Xi_IEt_VB_pow_decay_loglog_", length(X_vars), ".RData")
#-----------------------------------------
#-----------------------------------------
set.seed(1234)
#-----------------------------------------
# Load data
#-----------------------------------------
load("./data/Kenya_Score_Data_r.RData")
# Ken.Site <- Site
# Ken.G    <- G.mat

Kenya_Score_Data <- Kenya_Score_Data[Kenya_Score_Data$Village_ID %nin% c("KEN212", "KEN199", "KEN086"), ]
ind        <- which(Site$Village_ID %in% c("KEN212", "KEN199", "KEN086"))
Ken.Site       <- Site[Site$Village_ID %nin% c("KEN212", "KEN199", "KEN086"), ]
Kenya.Dist.c <- Kenya.Dist.c[-ind, -ind]
Ken.G        <- G.mat[-ind, -ind]


load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Tan.G    <- G.mat
Site     <- rbind(Ken.Site, Tan.Site)
# load("./data/sim.Cov.Data.RData")
#-----------------------------------------
# Define indices for regions
#-----------------------------------------
region_flags <- c("Northwestern", "Northeastern", "Southern", "Western", "Eastern")
Ken_indices  <- lapply(region_flags[1:3], function(f) which(Ken.Site$flag == f))
Tan_indices  <- lapply(region_flags[4:5], function(f) which(Tan.Site$flag == f))
Score_Data   <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
#-----------------------------------------
Y_vars <- c("Prevalence")
YEAR   <- unique(Score_Data$Year) %>% sort()
date.time <- data.frame(time.index = 1:length(YEAR),
                        time.scale = seq(0, 1, , length(YEAR)),
                        Year       = YEAR)
Score_Data$true.Prevalence <- Score_Data$Prevalence
Score_Data      <- Score_Data %>% left_join(date.time, by = c("Year"))
setDT(Site)
#-----------------------------------------
# Beta estimation per region
#-----------------------------------------
Beta.list    <- vector("list", length(region_flags))
G.basic.data <- vector("list", length(region_flags))
for(r in 1:length(region_flags)) {
  if(r <=3){
    Temp <- Ken.Site[Ken_indices[[r]], c("Village_ID","LAT","LON")] %>%
      left_join(Score_Data, by = c("Village_ID"))
  }else{
    Temp <- Tan.Site[Tan_indices[[r - 3]], c("Village_ID","LAT","LON")] %>%
      left_join(Score_Data, by = c("Village_ID"))
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
    data = subset(Score_Data, flag == f),
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
    obs.sigma.sq = list(mu.sigma.sq = 1e1, a = 2, b = 1),
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
                      CV                = FALSE,
                      Object            = "Flag",
                      transf.Response   = "loglog",
                      plot              = TRUE,
                      Ne                = 300,
                      tol.real          = 1e-5,
                      itMin             = 2e2,
                      itMax             = 1e3,
                      Obj.Seq           = 1)
end.time <- proc.time()
print(end.time - start.time)
run.time <- start.time - end.time
save(CV_Ranalysis, run.time, file = Tab.Name)











