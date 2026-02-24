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
Competing.method <- readxl::read_xlsx(paste0(root.table, "/misspecified_x_competing_methods_250.xlsx"))
save.Tab.from1   <- "./result/Simulation_300/misspecified_x_random_JSTVC_n_250_Ne_300/"
save.Tab.from2   <- "./result/Simulation_300/misspecified_x_random_JSTVC_IE_n_250_Ne_300/"
save.Tab.to <- paste0(root.table, "/misspecified_x_all_250.xlsx")


Competing.method <- readxl::read_xlsx(paste0(root.table, "/misspecified_x_competing_methods_298.xlsx"))
save.Tab.from1   <- "./result/Simulation_300/misspecified_x_random_JSTVC_n_250_Ne_300/"
save.Tab.from2   <- "./result/Simulation_300/misspecified_x_random_JSTVC_IE_n_250_Ne_300/"
save.Tab.to <- paste0(root.table, "/misspecified_x_all_298.xlsx")

#-----------------------------------------
#-----------------------------------------
beta_true <- c(5, rep(-1, 10))
JSTVC_beta_mean <- JSTVC_beta_var <- JSTVC_IE_beta_mean <- JSTVC_IE_beta_var <- NULL
# JSTVC
JSTVC_files <- list.files(save.Tab.from1)
JSTVC_files <- JSTVC_files[order(as.numeric(sub("sim_(\\d+)\\.RData", "\\1", JSTVC_files)))]
for(iter in 1:length(JSTVC_files)){
  load(paste0(save.Tab.from1, JSTVC_files[iter]))
  JSTVC_beta_mean <- rbind(JSTVC_beta_mean, t(beta))
  JSTVC_beta_var <- rbind(JSTVC_beta_var, t(sigma.sq))
}

# JSTVC_IE
JSTVC_IE_files <- list.files(save.Tab.from2)
JSTVC_IE_files <- JSTVC_IE_files[order(as.numeric(sub("sim_(\\d+)\\.RData", "\\1", JSTVC_IE_files)))]
for(iter in 1:(length(JSTVC_IE_files))){
  load(paste0(save.Tab.from2, JSTVC_IE_files[iter]))

  JSTVC_IE_beta_mean <- rbind(JSTVC_IE_beta_mean, t(beta))
  JSTVC_IE_beta_var <- rbind(JSTVC_IE_beta_var, t(sigma.sq))
}

JSTVC     <- calc_performance(JSTVC_beta_mean, beta_true)
JSTVC_ie  <- calc_performance(JSTVC_IE_beta_mean, beta_true[1:9])

JSTVC$Details$JSTVC.EC   <-  round(colMeans(compute_coverage(JSTVC_beta_mean, JSTVC_beta_var, beta_true))*100, 5)
JSTVC_ie$Details$JSTVC.EC   <-  round(colMeans(compute_coverage(JSTVC_IE_beta_mean, JSTVC_IE_beta_var, beta_true[1:9]))*100, 5)





create_wide_table <- function(Competing.method,
                              JSTVC_ie = NULL,
                              JSTVC,
                              metric = "Bias",
                              models = c("JSTVC-ξ", "Sub-AR1", "Arm-AR1", "JSTVC_IE", "JSTVC"),
                              n_beta = 11,
                              digits = 3,
                              Opt    = TRUE) {
  # Generate Beta column names
  beta_names <- paste0("Beta", 1:n_beta)

  # Initialize the data.frame with Metric and Model columns
  da_wide <- data.frame(
    Metric = metric,
    Model = models
  )
  cols <- grep(paste0("\\.", tolower(metric), "$"), names(Competing.method), value = TRUE)
  beta_matrix <- NULL
  for(l in 1:length(cols)){
    x <- Competing.method[[ cols[l] ]]
    beta_matrix <- rbind(beta_matrix, t(formatC(x, format = "f", digits = digits)))
    # beta_matrix <- rbind(beta_matrix, t(Competing.method[, cols[l]]))
  }
  # Now add JSTVC$Details

  if(!is.null(JSTVC_ie)){
    cols_jstvc <- grep((metric), names(JSTVC_ie$Details), value = TRUE) #tolower
    jstvc_ie_matrix <- cbind(do.call(rbind, lapply(cols_jstvc, function(col) {
      t(formatC(JSTVC_ie$Details[, col], format = "f", digits = digits))
    })), "", "")
    beta_matrix <- rbind(beta_matrix, jstvc_ie_matrix)
  }
  cols_jstvc <- grep((metric), names(JSTVC$Details), value = TRUE) #tolower
  jstvc_matrix <- do.call(rbind, lapply(cols_jstvc, function(col) {
    t(formatC(JSTVC$Details[, col], format = "f", digits = digits))
  }))
  beta_matrix <- rbind(beta_matrix,  jstvc_matrix)


  # Assign column names
  colnames(beta_matrix) <- beta_names

  # Add Beta columns to the data.frame
  da_wide[, beta_names] <- beta_matrix

  if(Opt){
    for (col in 3:(n_beta + 2)) {
      numeric_vals <- abs(as.numeric(da_wide[[col]]))  # ensure values are numeric
      max_val <- min(abs(numeric_vals), na.rm = TRUE)
      da_wide[[col]] <- ifelse(abs(numeric_vals) == max_val,
                               sprintf("\\textbf{%.4f}", numeric_vals),
                               sprintf("%.4f", numeric_vals))
    }
  }
  return(da_wide)
}

models <- c("$\\mathrm{JSTVC}_{-\\xi}$",
            "Sub-$\\mathrm{AR}_{1}$",
            "Arm-$\\mathrm{AR}_{1}$",
            "$\\mathrm{JSTVC}_{-\\mathrm{IE}}$",
            "JSTVC")

Bias <- create_wide_table(Competing.method, JSTVC_ie = JSTVC_ie, JSTVC, metric = "Bias", models = models, digits = 5)
SD   <- create_wide_table(Competing.method, JSTVC_ie= JSTVC_ie, JSTVC, metric = "SD", models = models, digits = 5)
MSE  <- create_wide_table(Competing.method, JSTVC_ie= JSTVC_ie, JSTVC, metric = "MSE", models = models, digits = 5)
EC   <- create_wide_table(Competing.method, JSTVC_ie= JSTVC_ie, JSTVC, metric = "EC",  models = models, digits = 2, Opt = F)



# writexl::write_xlsx(EC, path = "./Result/EC2.xlsx")

writexl::write_xlsx(rbind(Bias, SD, MSE), path = save.Tab.to)
# Table1
# Bias
# SD
# MSE

