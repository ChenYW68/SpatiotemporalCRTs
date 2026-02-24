rm(list = ls())
source("./LoadPackages/RDependPackages.R")
set.seed(12345)
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site

load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Site     <- rbind(Ken.Site, Tan.Site)
Score_Data <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
#3. others
X_vars <-c(c("Intercept")
           , paste0("CWT_", 1:4)
           , paste0("SBT_", 1:4)
           , "IEt.CWT"
           , "IEt.SBT")
JSTVC_xi <- paste0("Prevalence ~ -1 +  ", paste0(X_vars, collapse = " + "))

Sub_AR1 <- paste0("Prevalence ~ -1 + f(Year, model = 'ar1') + f(flag, model = 'iid') + ",
                  paste0(X_vars, collapse = " + "))

Arm_AR1 <- paste0("Prevalence ~ -1 + f(Year, model = 'ar1') + f(Study_Arm, model = 'iid') + ",
                  paste0(X_vars, collapse = " + "))
#3.1 JSTVC_xi
JSTVC_xi.fit  <- inla(as.formula(JSTVC_xi),
                      family = "gaussian",
                      data   = Score_Data,
                      verbose = FALSE,
                      control.predictor = list(compute = TRUE))

JSTVC_xi.beta <- JSTVC_xi.fit$summary.fixed$mean[-1]
#3.2 Sub_AR1
Sub_AR1.fit  <- inla(as.formula(Sub_AR1),
                     family = "gaussian",
                     data   = Score_Data,
                     verbose = FALSE,
                     control.predictor = list(compute = TRUE))

Sub_AR1.beta <- Sub_AR1.fit$summary.fixed$mean[-1]

#3.3 Arm_AR1
Arm_AR1.fit  <- inla(as.formula(Arm_AR1),
                     family = "gaussian",
                     data   = Score_Data,
                     verbose = FALSE,
                     control.predictor = list(compute = TRUE))

Arm_AR1.beta <- Arm_AR1.fit$summary.fixed$mean[-1]

save(JSTVC_xi.beta, Sub_AR1.beta, Arm_AR1.beta, file = "./result/case/three_competing_models.RData")



