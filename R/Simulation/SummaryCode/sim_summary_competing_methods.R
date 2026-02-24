rm(list = ls())
#-----------------------------------------
source("./LoadPackages/RDependPackages.R")
source(normalizePath("./JSTVC/R/util.R"))
#-----------------------------------------
root.table <- "./result/summary"
if (!dir.exists(root.table)) {
  dir.create(root.table, recursive = TRUE)
}
#-----------------------------------------
#-----------------------------------------
save.Tab.from  <- "./result/Simulation_300/misspecified_x_random_competing_n_250/"
save.Tab.to <- paste0(root.table, "/misspecified_x_competing_methods_250.xlsx")

# save.Tab.from  <- "./result/Simulation_300/misspecified_x_random_competing_n_298/"
# save.Tab.to <- paste0(root.table, "/misspecified_x_competing_methods_298.xlsx")
#-----------------------------------------
#-----------------------------------------
beta_true <- c(5, rep(-1, 10))
nc <- 300;
nonEffect.beta_mean <- subRegion.beta_mean <- subArm.beta_mean <- NULL
nonEffect.beta_var <- subRegion.beta_var <- subArm.beta_var<- NULL
file <- list.files(save.Tab.from)
file <- file[order(as.numeric(sub("sim_(\\d+)\\.RData", "\\1", file)))]
# file <- file#[1:nc]
for(iter in 1:(length(file))){
  load(paste0(save.Tab.from, file[iter]))
  nonEffect.beta_mean <- rbind(nonEffect.beta_mean, t(JSTVC.nonEffect.beta))
  subRegion.beta_mean <- rbind(subRegion.beta_mean, t(JSTVC.subRegion.beta))
  subArm.beta_mean <- rbind(subArm.beta_mean, t(JSTVC.subArm.beta))

  nonEffect.beta_var <- rbind(nonEffect.beta_var, t(JSTVC.nonEffect.sd^2))
  subRegion.beta_var <- rbind(subRegion.beta_var, t(JSTVC.subRegion.sd^2))
  subArm.beta_var <- rbind(subArm.beta_var, t(JSTVC.subArm.sd^2))
}


nonEffect  <- calc_performance(nonEffect.beta_mean, beta_true)
subRegion  <- calc_performance(subRegion.beta_mean, beta_true)
subArm  <- calc_performance(subArm.beta_mean, beta_true)

nonEffect.coverage      <- compute_coverage(nonEffect.beta_mean, nonEffect.beta_var, beta_true)
subRegion.coverage      <- compute_coverage(subRegion.beta_mean, subRegion.beta_var, beta_true)
subArm.coverage      <- compute_coverage(subArm.beta_mean, subArm.beta_var, beta_true)

nonEffect.coverage_prob   <- round(colMeans(nonEffect.coverage)*100, 5)
subRegion.coverage_prob   <- round(colMeans(subRegion.coverage)*100, 5)
subArm.coverage_prob     <- round(colMeans(subArm.coverage)*100, 5)


da <- data.frame(beta = nonEffect$Details$Coefficient,
                 nonEffect.bias = nonEffect$Details$Bias,
                 subRegion.bias = subRegion$Details$Bias,
                 subArm.bias = subArm$Details$Bias,

                 nonEffect.sd = nonEffect$Details$SD,
                 subRegion.sd = subRegion$Details$SD,
                 subArm.sd = subArm$Details$SD,

                 nonEffect.mse = nonEffect$Details$MSE,
                 subRegion.mse = subRegion$Details$MSE,
                 subArm.mse    = subArm$Details$MSE,

                 nonEffect.ec = nonEffect.coverage_prob,
                 subRegion.ec = subRegion.coverage_prob,
                 subArm.ec = subArm.coverage_prob
)
writexl::write_xlsx(da, path = save.Tab.to)

# da
# nonEffect.coverage_prob
# length(file)
