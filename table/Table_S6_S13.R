rm(list = ls())
source("./LoadPackages/RDependPackages.R")
source(normalizePath("./JSTVC/R/util.R"))
beta_true <- c(5, rep(-1, 10))

#Output source
#-----------------------------------------
root.table <- "./result/summary"
if (!dir.exists(root.table)) {
  dir.create(root.table, recursive = TRUE)
}
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#--Results from the Scenario Using the Gneiting Spatiotemporal Covariance Structure
#1. Competing models
write_competing_summary(
  save.Tab.from = "./result/Simulation/random_competing_n_250/",
  save.Tab.to = file.path(root.table, "random_competing_methods_250.xlsx"),
  beta_true = beta_true
)
write_competing_summary(
  save.Tab.from = "./result/Simulation/random_competing_n_298/",
  save.Tab.to = file.path(root.table, "random_competing_methods_298.xlsx"),
  beta_true = beta_true
)

#2. Merging with those from JSTVC
#Table S6, and EC of S8
summarize_one_setting(
  competing_path = paste0(root.table, "/random_competing_methods_250.xlsx"),
  save.Tab.from1 = "./result/Simulation/random_JSTVC_n_250/",
  save.Tab.from2 = "./result/Simulation/random_JSTVC_IE_n_250/",
  save.Tab.to    = paste0(root.table, "/Table_S6_S8_random_all_250.xlsx")
)
#Table S7, and EC of S8
summarize_one_setting(
  competing_path = paste0(root.table, "/random_competing_methods_298.xlsx"),
  save.Tab.from1 = "./result/Simulation/random_JSTVC_n_298/",
  save.Tab.from2 = "./result/Simulation/random_JSTVC_IE_n_298/",
  save.Tab.to    = paste0(root.table, "/Table_S7_S8_random_all_298.xlsx")
)





#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#--Results from the Scenario Using Low-Rank Basis Expansion (Smoothed)
#1. Competing models
write_competing_summary(
  save.Tab.from = "./result/Simulation/smoothed_competing_n_250/",
  save.Tab.to = file.path(root.table, "smoothed_competing_methods_250.xlsx"),
  beta_true = beta_true
)
write_competing_summary(
  save.Tab.from = "./result/Simulation/smoothed_competing_n_298/",
  save.Tab.to = file.path(root.table, "smoothed_competing_methods_298.xlsx"),
  beta_true = beta_true
)


#2. Merging with those from JSTVC
#Table S9, and EC of S11
summarize_one_setting(
  competing_path = paste0(root.table, "/smoothed_competing_methods_250.xlsx"),
  save.Tab.from1 = "./result/Simulation/smoothed_JSTVC_n_250/",
  save.Tab.from2 = "./result/Simulation/smoothed_JSTVC_IE_n_250/",
  save.Tab.to    = paste0(root.table, "/Table_S9_S11_smoothed_all_250.xlsx")
)
#Table S10, and EC of S11
summarize_one_setting(
  competing_path = paste0(root.table, "/smoothed_competing_methods_298.xlsx"),
  save.Tab.from1 = "./result/Simulation/smoothed_JSTVC_n_298/",
  save.Tab.from2 = "./result/Simulation/smoothed_JSTVC_IE_n_298/",
  save.Tab.to    = paste0(root.table, "/Table_S10_S11_smoothed_all_298.xlsx")
)



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#--Results from the Scenario for Misspecification of the Decay Function
#Competing models
write_competing_summary(
  save.Tab.from = "./result/Simulation/misspecified_x_random_competing_n_250/",
  save.Tab.to = file.path(root.table, "misspecified_x_random_competing_methods_250.xlsx"),
  beta_true = beta_true
)
write_competing_summary(
  save.Tab.from = "./result/Simulation/misspecified_x_random_competing_n_298/",
  save.Tab.to = file.path(root.table, "misspecified_x_random_competing_methods_298.xlsx"),
  beta_true = beta_true
)

#Table S12
summarize_one_setting(
  competing_path = paste0(root.table, "/misspecified_x_random_competing_methods_250.xlsx"),
  save.Tab.from1 = "./result/Simulation/misspecified_x_random_JSTVC_n_250/",
  save.Tab.from2 = "./result/Simulation/misspecified_x_random_JSTVC_IE_n_250/",
  save.Tab.to    = paste0(root.table, "/Table_S12_misspecified_x_all_250.xlsx")
)
#Table S13
summarize_one_setting(
  competing_path = paste0(root.table, "/misspecified_x_random_competing_methods_298.xlsx"),
  save.Tab.from1 = "./result/Simulation/misspecified_x_random_JSTVC_n_298/",
  save.Tab.from2 = "./result/Simulation/misspecified_x_random_JSTVC_IE_n_298/",
  save.Tab.to    = paste0(root.table, "/Table_S13_misspecified_x_all_298.xlsx")
)

