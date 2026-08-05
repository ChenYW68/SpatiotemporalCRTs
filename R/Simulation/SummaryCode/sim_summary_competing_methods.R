rm(list = ls())
#-----------------------------------------
source("./LoadPackages/RDependPackages.R")
Rcpp::sourceCpp("./JSTVC/src/util_c.cpp")
source(normalizePath("./JSTVC/R/util.R"))
#-----------------------------------------
root.table <- "./result/summary"
if (!dir.exists(root.table)) {
  dir.create(root.table, recursive = TRUE)
}
#-----------------------------------------
beta_true <- c(5, rep(-1, 10))



write_competing_summary(
  save.Tab.from = "./result/Simulation_300/misspecified_x_random_competing_n_250/",
  save.Tab.to = file.path(root.table, "misspecified_x_competing_methods_250.xlsx"),
  beta_true = beta_true
)

# write_competing_summary(
#   save.Tab.from = "./result/Simulation_300/misspecified_x_random_competing_n_298/",
#   save.Tab.to = file.path(root.table, "misspecified_x_competing_methods_298.xlsx"),
#   beta_true = beta_true
# )
