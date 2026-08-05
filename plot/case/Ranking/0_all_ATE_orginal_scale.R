rm(list = ls())

source("./LoadPackages/RDependPackages.R")
source(normalizePath("./JSTVC/R/util.R"))

set.seed(12345)

posterior_sample_size <- 5000

load("./data/Kenya_Score_Data_r.RData")
Kenya.Score.Data <- as.data.frame(Kenya_Score_Data)

load("./data/Tanzania_Score_Data_r.RData")
Tanzania.Score.Data <- as.data.frame(Tanzania_Score_Data)

Score_Data <- rbind(Kenya.Score.Data, Tanzania.Score.Data)

output_dir <- "./result/case/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

standard_treatment_vars <- c(
  "Intercept",
  paste0("CWT_", 1:4),
  paste0("SBT_", 1:4),
  "IEt.CWT",
  "IEt.SBT"
)
no_ie_treatment_vars <- c("Intercept", paste0("CWT_", 1:4), paste0("SBT_", 1:4))
exp_treatment_vars <- c(
  "Intercept",
  paste0("CWT_", 1:4),
  paste0("SBT_", 1:4),
  "IEt.CWT.exp",
  "IEt.SBT.exp"
)
non_decay_treatment_vars <- c(
  "Intercept",
  paste0("CWT_", 1:4),
  paste0("SBT_", 1:4),
  paste0("sCWT_", 1:2),
  paste0("sSBT_", 1:2)
)
mixed_treatment_vars <- c("Intercept", "CWT", "SBT")
spie_treatment_vars <- function(neighbour_count) {
  c(
    "Intercept",
    paste0("CWT_", 1:4),
    paste0("SBT_", 1:4),
    "IEt.CWT",
    "IEt.SBT",
    paste0("IEs.CWT.sp.Neigh.500.", neighbour_count),
    paste0("IEs.SBT.sp.Neigh.500.", neighbour_count)
  )
}

inverse_cloglog <- function(eta) {
  1 - exp(-exp(eta))
}

cloglog_transform <- function(x) {
  log(-log(1 - x))
}

identity_transform <- function(eta) {
  eta
}

compute_true_prevalence <- function(score_data) {
  score_data |>
    dplyr::select(Village_ID, Year, Prevalence, Study_Arm) |>
    dplyr::group_by(Year, Study_Arm) |>
    dplyr::summarise(Prevalence = mean(Prevalence, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = Study_Arm,
      values_from = Prevalence
    ) |>
    as.data.frame() |>
    (\(x) x[, -1, drop = FALSE])()
}

compute_spatiotemporal_ate <- function(score_data,
                                       fit_object,
                                       treatment_names,
                                       transform_fn,
                                       beta_sampler = c("vb", "mcmc")) {
  beta_sampler <- match.arg(beta_sampler)
  arm_profiles <- extract_arm_profiles(score_data, treatment_names)
  year_values <- sort(unique(arm_profiles$Year))
  arm_ids <- sort(unique(arm_profiles$Study_Arm))
  reference_profile <- make_reference_regime(year_values, treatment_names)
  random_effect_draw_indices <- NULL
  random_effect_draws_override <- NULL
  obs_sigma_sq_draws <- NULL

  beta_draws <- switch(
    beta_sampler,
    vb = {
      random_effect_draw_indices <- sample_random_effect_draw_indices(
        fit_object,
        n_draws = posterior_sample_size
      )
      obs_sigma_sq_draws <- sample_obs_sigma_sq_draws_vb(
        fit_object,
        n_draws = posterior_sample_size
      )
      sample_beta_draws_vb(
        fit_object,
        n_draws = posterior_sample_size
      )
    },
    mcmc = {
      random_effect_draws_override <- extract_random_effect_draws_mcmc(
        fit_object,
        n_draws = posterior_sample_size
      )
      obs_sigma_sq_draws <- sample_obs_sigma_sq_draws_mcmc(
        fit_object,
        n_draws = posterior_sample_size
      )
      sample_beta_draws_mcmc(
        fit_object,
        n_draws = posterior_sample_size
      )
    }
  )

  reference_site_outcomes <- compute_regime_site_outcomes(
    regime_profile = reference_profile,
    fit_object = fit_object,
    beta_draws = beta_draws,
    obs_sigma_sq_draws = obs_sigma_sq_draws,
    treatment_names = treatment_names,
    transform_fn = transform_fn,
    random_effect_draw_indices = random_effect_draw_indices,
    random_effect_draws_override = random_effect_draws_override
  )
  reference_mean <- average_site_outcomes(reference_site_outcomes)

  ate_draws <- array(
    0,
    dim = c(length(year_values), length(arm_ids), nrow(beta_draws)),
    dimnames = list(
      as.character(year_values),
      paste0("Arm_", arm_ids),
      paste0("Draw_", seq_len(nrow(beta_draws)))
    )
  )

  for (arm_index in seq_along(arm_ids)) {
    arm_profile <- arm_profiles[arm_profiles$Study_Arm == arm_ids[arm_index], c("Year", treatment_names), drop = FALSE]
    arm_site_outcomes <- compute_regime_site_outcomes(
      regime_profile = arm_profile,
      fit_object = fit_object,
      beta_draws = beta_draws,
      obs_sigma_sq_draws = obs_sigma_sq_draws,
      treatment_names = treatment_names,
      transform_fn = transform_fn,
      random_effect_draw_indices = random_effect_draw_indices,
      random_effect_draws_override = random_effect_draws_override
    )
    arm_mean <- average_site_outcomes(arm_site_outcomes)
    ate_draws[, arm_index, ] <- arm_mean - reference_mean
  }

  list(
    arm_profiles = arm_profiles,
    reference_profile = reference_profile,
    ATE_draws = ate_draws,
    ATE = apply(ate_draws, c(1, 2), mean)
  )
}

fit_competing_inla_models <- function(score_data,
                                      transform_fn = cloglog_transform) {
  score_data <- as.data.frame(score_data)
  score_data$Prevalence <- transform_fn(score_data$Prevalence)

  x_vars <- c(
    "Intercept",
    paste0("CWT_", 1:4),
    paste0("SBT_", 1:4),
    "IEt.CWT",
    "IEt.SBT"
  )

  jstvc_xi_formula <- paste0("Prevalence ~ -1 + ", paste0(x_vars, collapse = " + "))
  sub_ar1_formula <- paste0(
    "Prevalence ~ -1 + f(Year, model = 'ar1') + f(flag, model = 'iid') + ",
    paste0(x_vars, collapse = " + ")
  )
  arm_ar1_formula <- paste0(
    "Prevalence ~ -1 + f(Year, model = 'ar1') + f(Study_Arm, model = 'iid') + ",
    paste0(x_vars, collapse = " + ")
  )

  list(
    JSTVC_xi.fit = inla(
      as.formula(jstvc_xi_formula),
      family = "gaussian",
      data = score_data,
      verbose = FALSE,
      control.predictor = list(compute = TRUE)
    ),
    Sub_AR1.fit = inla(
      as.formula(sub_ar1_formula),
      family = "gaussian",
      data = score_data,
      verbose = FALSE,
      control.predictor = list(compute = TRUE)
    ),
    Arm_AR1.fit = inla(
      as.formula(arm_ar1_formula),
      family = "gaussian",
      data = score_data,
      verbose = FALSE,
      control.predictor = list(compute = TRUE)
    )
  )
}

compute_inla_identity_ate <- function(score_data,
                                      fit_object,
                                      treatment_names,
                                      transform_fn = inverse_cloglog) {
  extract_random_effect_mean <- function(effect_name, ids) {
    if (!effect_name %in% names(fit_object$summary.random)) {
      return(rep(0, length(ids)))
    }

    random_summary <- fit_object$summary.random[[effect_name]]
    random_ids <- random_summary[[1]]
    random_mean <- random_summary$mean
    matched_effects <- random_mean[match(ids, random_ids)]
    matched_effects[is.na(matched_effects)] <- 0
    matched_effects
  }

  arm_profiles <- extract_arm_profiles(score_data, treatment_names)
  year_values <- sort(unique(arm_profiles$Year))
  arm_ids <- sort(unique(arm_profiles$Study_Arm))
  reference_profile <- make_reference_regime(year_values, treatment_names)

  fixed_effects <- stats::setNames(
    fit_object$summary.fixed$mean,
    rownames(fit_object$summary.fixed)
  )[treatment_names]
  current_design <- as.matrix(score_data[, treatment_names, drop = FALSE])
  current_treatment_effect <- as.vector(current_design %*% fixed_effects)
  observed_arm_effect <- extract_random_effect_mean("Study_Arm", score_data$Study_Arm)
  baseline_eta <- fit_object$summary.fitted.values$mean - current_treatment_effect - observed_arm_effect

  ate_matrix <- matrix(
    0,
    nrow = length(year_values),
    ncol = length(arm_ids),
    dimnames = list(as.character(year_values), paste0("Arm_", arm_ids))
  )

  for (arm_index in seq_along(arm_ids)) {
    arm_id <- arm_ids[arm_index]
    arm_profile <- arm_profiles[arm_profiles$Study_Arm == arm_id, c("Year", treatment_names), drop = FALSE]

    for (year_index in seq_along(year_values)) {
      year_value <- year_values[year_index]
      row_index <- score_data$Year == year_value

      arm_x <- as.numeric(arm_profile[arm_profile$Year == year_value, treatment_names, drop = FALSE][1, ])
      ref_x <- as.numeric(reference_profile[reference_profile$Year == year_value, treatment_names, drop = FALSE][1, ])

      eta_arm <- baseline_eta[row_index] + as.vector(crossprod(arm_x, fixed_effects))
      eta_ref <- baseline_eta[row_index] + as.vector(crossprod(ref_x, fixed_effects))
      ate_matrix[year_index, arm_index] <- mean(
        transform_fn(eta_arm) - transform_fn(eta_ref),
        na.rm = TRUE
      )
    }
  }

  list(
    arm_profiles = arm_profiles,
    reference_profile = reference_profile,
    ATE = ate_matrix
  )
}

compute_country_ate <- function(score_data,
                                fit_file,
                                treatment_names = standard_treatment_vars,
                                transform_fn = inverse_cloglog,
                                beta_sampler = "vb") {
  if (!file.exists(fit_file)) {
    stop(
      paste0(
        "Loaded data does not exist locally: ",
        fit_file,
        ". Please try again after running the corresponding model."
      )
    )
  }
  load(fit_file)
  compute_spatiotemporal_ate(
    score_data = score_data,
    fit_object = CV_Ranalysis,
    treatment_names = treatment_names,
    transform_fn = transform_fn,
    beta_sampler = beta_sampler
  )
}

# extract_vb_beta_mean <- function(fit_file) {
#   load(fit_file)
#   as.vector(CV_Ranalysis$update.Para.List[[1]]$beta$mu.beta[-1, 1])
# }

True. <- compute_true_prevalence(Score_Data)
True.Kenya <- compute_true_prevalence(Kenya.Score.Data)
True.Tanzania <- compute_true_prevalence(Tanzania.Score.Data)

#1. JSTVC_IE
print("Evaluating ATEs from JSTVC without IEs ...")
JSTVC_IE.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_VB_pow_decay_loglog_9.RData",
  treatment_names = no_ie_treatment_vars,
  transform_fn = inverse_cloglog
)
JSTVC_IE.ATE <- JSTVC_IE.results$ATE

#2.1 JSTVC
print("Evaluating ATEs using JSTVC with standard configurations ...")
JSTVC.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData",
  transform_fn = inverse_cloglog
)
JSTVC.ATE <- JSTVC.results$ATE

#2.2 JSTVC exp
print("Evaluating ATEs using JSTVC with the exponential decay function ...")
JSTVC_exp.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_VB_pow_decay_exp_11.RData",
  treatment_names = exp_treatment_vars,
  transform_fn = inverse_cloglog
)
JSTVC_exp.ATE <- JSTVC_exp.results$ATE

#2.3 JSTVC with original scale
print("Evaluating ATEs using JSTVC with the untransformed data ...")
JSTVC_or.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_VB_pow_decay_Untransform_11.RData",
  transform_fn = identity_transform
)
JSTVC_or.ATE <- JSTVC_or.results$ATE

#2.4 JSTVC with Kenya
print("Evaluating ATEs using JSTVC with the Kenya's data only ...")
JSTVC_Kenya.results <- compute_country_ate(
  score_data = Kenya.Score.Data,
  fit_file = "./result/case/Kenya_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData",
  transform_fn = inverse_cloglog
)
JSTVC_Kenya.ATE <- JSTVC_Kenya.results$ATE

#2.5 JSTVC with Tanzania
print("Evaluating ATEs using JSTVC with the Tanzania's data only ...")
JSTVC_Tanzania.results <- compute_country_ate(
  score_data = Tanzania.Score.Data,
  fit_file = "./result/case/Tanzania_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData",
  transform_fn = inverse_cloglog
)
JSTVC_Tanzania.ATE <- JSTVC_Tanzania.results$ATE


#2.6 JSTVC with MCMC
print("Evaluating ATEs using JSTVC with the MCMC implementation ...")
JSTVC_mcmc.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_MCMC_pow_decay_loglog_11.RData",
  beta_sampler = "mcmc",
  transform_fn = inverse_cloglog
)
JSTVC_mcmc.ATE <- JSTVC_mcmc.results$ATE

#2.7 JSTVC without outliers
print("Evaluating ATEs using JSTVC without outliers ...")
JSTVC_outlier.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/Delete_Outliers_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData",
  transform_fn = inverse_cloglog
)
JSTVC_outlier.ATE <- JSTVC_outlier.results$ATE

#2.8 JSTVC with mixed DE and IE
print("Evaluating ATEs using JSTVC with mixed DEs and IEs ...")
JSTVC_mixed.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_Mixed_VB_pow_decay_loglog_3.RData",
  treatment_names = mixed_treatment_vars,
  transform_fn = inverse_cloglog
)
JSTVC_mixed.ATE <- JSTVC_mixed.results$ATE

#2.9.1 JSTVC with spatial IEs (range = 10km)
print("Evaluating ATEs using JSTVC with spatiotemporal IEs (range = 10km) ...")
JSTVC_spIE_10.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_IEs_VB_pow_decay_loglog_13_10.RData",
  treatment_names = spie_treatment_vars(10),
  transform_fn = inverse_cloglog
)
JSTVC_spIE_10.ATE <- JSTVC_spIE_10.results$ATE

#2.9.2 JSTVC with spatial IEs (range = 30km)
print("Evaluating ATEs using JSTVC with spatiotemporal IEs (range = 30km) ...")
JSTVC_spIE_30.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_IEs_VB_pow_decay_loglog_13_30.RData",
  treatment_names = spie_treatment_vars(30),
  transform_fn = inverse_cloglog
)
JSTVC_spIE_30.ATE <- JSTVC_spIE_30.results$ATE

#2.9.3 JSTVC with spatial IEs (range = 50km)
print("Evaluating ATEs using JSTVC with spatiotemporal IEs (range = 50km) ...")
JSTVC_spIE_50.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_IEs_VB_pow_decay_loglog_13_50.RData",
  treatment_names = spie_treatment_vars(50),
  transform_fn = inverse_cloglog
)
JSTVC_spIE_50.ATE <- JSTVC_spIE_50.results$ATE

#2.10 JSTVC with non-decay_IEs
print("Evaluating ATEs using JSTVC with time-varying coefficients for IEs ...")
JSTVC_non_decay.results <- compute_country_ate(
  score_data = Score_Data,
  fit_file = "./result/case/JSTVC_Xi_IEt_VB_non_decay_loglog_13.RData",
  treatment_names = non_decay_treatment_vars,
  transform_fn = inverse_cloglog
)
JSTVC_non_decay.ATE <- JSTVC_non_decay.results$ATE


#3 Competing methods
competing_fits <- fit_competing_inla_models(
  Score_Data,
  transform_fn = cloglog_transform
)

print("Evaluating ATEs using JSTVC without xi_t(s) from INLA ...")
JSTVC_xi.results <- compute_inla_identity_ate(
  score_data = Score_Data,
  fit_object = competing_fits$JSTVC_xi.fit,
  treatment_names = standard_treatment_vars,
  transform_fn = inverse_cloglog
)
JSTVC_xi.ATE <- JSTVC_xi.results$ATE

print("Evaluating ATEs using sub-AR1 from INLA ...")
Sub_AR1.results <- compute_inla_identity_ate(
  score_data = Score_Data,
  fit_object = competing_fits$Sub_AR1.fit,
  treatment_names = standard_treatment_vars,
  transform_fn = inverse_cloglog
)
Sub_AR1.ATE <- Sub_AR1.results$ATE

print("Evaluating ATEs using Arm-AR1 from INLA ...")
Arm_AR1.results <- compute_inla_identity_ate(
  score_data = Score_Data,
  fit_object = competing_fits$Arm_AR1.fit,
  treatment_names = standard_treatment_vars,
  transform_fn = inverse_cloglog
)
Arm_AR1.ATE <- Arm_AR1.results$ATE

ATE_by_model <- list(
  JSTVC_IE.ATE = JSTVC_IE.ATE,
  JSTVC.ATE = JSTVC.ATE,
  JSTVC_xi.ATE = JSTVC_xi.ATE,
  Sub_AR1.ATE = Sub_AR1.ATE,
  Arm_AR1.ATE = Arm_AR1.ATE,
  JSTVC_exp.ATE = JSTVC_exp.ATE,
  JSTVC_outlier.ATE = JSTVC_outlier.ATE,
  JSTVC_or.ATE = JSTVC_or.ATE,
  JSTVC_Kenya.ATE = JSTVC_Kenya.ATE,
  JSTVC_Tanzania.ATE = JSTVC_Tanzania.ATE,
  JSTVC_mcmc.ATE = JSTVC_mcmc.ATE,
  JSTVC_mixed.ATE = JSTVC_mixed.ATE,
  JSTVC_spIE_10.ATE = JSTVC_spIE_10.ATE,
  JSTVC_spIE_30.ATE = JSTVC_spIE_30.ATE,
  JSTVC_spIE_50.ATE = JSTVC_spIE_50.ATE,
  JSTVC_non_decay.ATE = JSTVC_non_decay.ATE
)

save(
  ATE_by_model,
  True.,
  True.Kenya,
  True.Tanzania,
  JSTVC_IE.results,
  JSTVC_IE.ATE,
  JSTVC.results,
  JSTVC.ATE,
  JSTVC_xi.results,
  JSTVC_xi.ATE,
  Sub_AR1.results,
  Sub_AR1.ATE,
  Arm_AR1.results,
  Arm_AR1.ATE,
  JSTVC_exp.results,
  JSTVC_exp.ATE,
  JSTVC_outlier.results,
  JSTVC_outlier.ATE,
  JSTVC_or.results,
  JSTVC_or.ATE,
  JSTVC_Kenya.results,
  JSTVC_Kenya.ATE,
  JSTVC_Tanzania.results,
  JSTVC_Tanzania.ATE,
  JSTVC_mcmc.results,
  JSTVC_mcmc.ATE,
  JSTVC_mixed.results,
  JSTVC_mixed.ATE,
  JSTVC_spIE_10.results,
  JSTVC_spIE_10.ATE,
  JSTVC_spIE_30.results,
  JSTVC_spIE_30.ATE,
  JSTVC_spIE_50.results,
  JSTVC_spIE_50.ATE,
  JSTVC_non_decay.results,
  JSTVC_non_decay.ATE,
  file = file.path(output_dir, "all_ATE_orginal_scale.RData")
)

print(lapply(ATE_by_model, colSums))
