rm(list = ls())

source("./LoadPackages/RDependPackages.R")

load("./data/Kenya_Score_Data_r.RData")
Kenya.Score.Data <- Kenya_Score_Data

load("./data/Tanzania_Score_Data_r.RData")
Tanzania.Score.Data <- Tanzania_Score_Data

Score_Data <- rbind(Kenya.Score.Data, Tanzania.Score.Data)
load("./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData")

output_dir <- "./result/case/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

treatment_vars <- c(
  "Intercept",
  paste0("CWT_", 1:4),
  paste0("SBT_", 1:4),
  "IEt.CWT",
  "IEt.SBT"
)
indirect_effect_vars <- c("IEt.CWT", "IEt.SBT")

transform_to_prevalence <- function(eta) {
  1 - exp(-exp(eta))
}

make_reference_regime <- function(years) {
  data.frame(
    Year = years,
    Intercept = 1,
    CWT_1 = 0,
    CWT_2 = 0,
    CWT_3 = 0,
    CWT_4 = 0,
    SBT_1 = 0,
    SBT_2 = 0,
    SBT_3 = 0,
    SBT_4 = 0,
    IEt.CWT = 0,
    IEt.SBT = 0
  )
}

extract_arm_profiles <- function(score_data, treatment_names) {
  score_data |>
    dplyr::select(Study_Arm, Year, dplyr::all_of(treatment_names)) |>
    dplyr::group_by(Study_Arm, Year) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(treatment_names), ~ .x[1]), .groups = "drop") |>
    dplyr::arrange(Study_Arm, Year)
}

make_de_only_profile <- function(arm_profile, ie_names) {
  de_profile <- arm_profile
  de_profile[, ie_names] <- 0
  de_profile
}

sample_beta_draws <- function(fit_object) {
  sigma <- fit_object$update.Para.List[[1]]$beta$sigma.sq
  if (length(dim(sigma)) != 2) {
    sigma <- as.matrix(sigma)
  }

  mvnfast::rmvn(
    n = dim(fit_object$fitted.Pred.ens[[1]])[3],
    mu = as.vector(fit_object$update.Para.List[[1]]$beta$mu.beta),
    sigma = sigma
  )
}

sample_obs_noise_draws <- function(fit_object, n_draws) {
  region_names <- names(fit_object$data)
  noise_draws <- vector("list", length(region_names))
  names(noise_draws) <- region_names

  for (region_index in seq_along(region_names)) {
    region_name <- region_names[region_index]
    n_years <- fit_object$data[[region_index]]$Nt
    n_sites <- fit_object$data[[region_index]]$n
    sigma_sq <- fit_object$update.Para.List[[region_name]]$obs.sigma.sq$mu.sigma.sq

    noise_draws[[region_name]] <- array(
      rnorm(n_years * n_sites * n_draws, mean = 0, sd = sqrt(sigma_sq)),
      dim = c(n_years, n_sites, n_draws),
      dimnames = list(
        rownames(fit_object$data[[region_index]]$Y_ts),
        colnames(fit_object$data[[region_index]]$Y_ts),
        paste0("Draw_", seq_len(n_draws))
      )
    )
  }

  noise_draws
}

compute_regime_site_outcomes <- function(regime_profile,
                                         fit_object,
                                         beta_draws,
                                         treatment_names,
                                         obs_noise_draws = NULL) {
  region_names <- names(fit_object$data)
  n_draws <- nrow(beta_draws)
  year_values <- sort(unique(regime_profile$Year))
  n_years <- length(year_values)

  site_outcomes <- vector("list", length(region_names))
  names(site_outcomes) <- region_names

  for (region_index in seq_along(region_names)) {
    region_name <- region_names[region_index]
    random_effect_draws <- fit_object$fitted.slope.Hv.Zg.ens[[paste0(region_name, ".Intercept")]][
      ,
      fit_object$data[[region_index]]$index.y,
      ,
      drop = FALSE
    ]
    n_sites <- dim(random_effect_draws)[2]
    site_ids <- colnames(fit_object$data[[region_index]]$Y_ts)
    if (is.null(site_ids)) {
      site_ids <- paste0(region_name, "_Site_", seq_len(n_sites))
    }

    site_outcomes[[region_name]] <- array(
      NA_real_,
      dim = c(n_years, n_sites, n_draws),
      dimnames = list(
        as.character(year_values),
        site_ids,
        paste0("Draw_", seq_len(n_draws))
      )
    )

    for (year_index in seq_len(n_years)) {
      year_value <- year_values[year_index]
      profile_row <- regime_profile[regime_profile$Year == year_value, treatment_names, drop = FALSE]

      if (nrow(profile_row) != 1) {
        stop("Each regime must provide exactly one treatment profile per year.")
      }

      x_vector <- as.numeric(profile_row[1, ])
      x_beta_draws <- as.vector(beta_draws %*% x_vector)
      eta_draws <- sweep(random_effect_draws[year_index, , , drop = FALSE][1, , ], 2, x_beta_draws, "+")
      if (!is.null(obs_noise_draws)) {
        eta_draws <- eta_draws + obs_noise_draws[[region_name]][year_index, , ]
      }
      # Potential outcomes are first mapped back to the original prevalence scale
      # at each site and only then averaged across sites.
      site_outcomes[[region_name]][year_index, , ] <- transform_to_prevalence(eta_draws)
    }
  }

  site_outcomes
}

average_site_outcomes <- function(site_outcomes) {
  year_values <- dimnames(site_outcomes[[1]])[[1]]
  draw_names <- dimnames(site_outcomes[[1]])[[3]]

  regime_mean <- matrix(
    0,
    nrow = length(year_values),
    ncol = length(draw_names),
    dimnames = list(year_values, draw_names)
  )

  total_sites <- sum(vapply(site_outcomes, function(x) dim(x)[2], numeric(1)))
  for (region_name in names(site_outcomes)) {
    regime_mean <- regime_mean + apply(site_outcomes[[region_name]], c(1, 3), sum) / total_sites
  }

  regime_mean
}

summarise_effect_draws <- function(effect_draws, arm_ids, year_values, effect_name) {
  effect_summary <- apply(effect_draws, c(1, 2), mean)
  effect_lower <- apply(effect_draws, c(1, 2), quantile, probs = 0.025)
  effect_upper <- apply(effect_draws, c(1, 2), quantile, probs = 0.975)
  effect_overall <- apply(effect_draws, c(2, 3), sum)
  effect_overall_summary <- data.frame(
    Study_Arm = arm_ids,
    Effect = apply(effect_overall, 1, mean),
    Lower = apply(effect_overall, 1, quantile, probs = 0.025),
    Upper = apply(effect_overall, 1, quantile, probs = 0.975)
  )
  names(effect_overall_summary)[2] <- effect_name

  effect_table <- expand.grid(
    Year = year_values,
    Study_Arm = arm_ids
  )
  effect_table$Effect <- as.vector(effect_summary)
  effect_table$Lower <- as.vector(effect_lower)
  effect_table$Upper <- as.vector(effect_upper)
  names(effect_table)[3] <- effect_name

  list(
    summary = effect_summary,
    lower = effect_lower,
    upper = effect_upper,
    overall = effect_overall,
    overall_summary = effect_overall_summary,
    table = effect_table
  )
}

build_decomposition_comparison <- function(one_step_draws,
                                           de_draws,
                                           ie_draws,
                                           arm_ids,
                                           year_values,
                                           effect_label) {
  de_plus_ie_draws <- de_draws + ie_draws
  comparison_diff_draws <- one_step_draws - de_plus_ie_draws

  comparison_table <- expand.grid(
    Year = year_values,
    Study_Arm = arm_ids
  )
  comparison_table[[paste0(effect_label, "_one_step")]] <- as.vector(apply(one_step_draws, c(1, 2), mean))
  comparison_table[[paste0(effect_label, "_DE_plus_IE")]] <- as.vector(apply(de_plus_ie_draws, c(1, 2), mean))
  comparison_table[[paste0(effect_label, "_diff")]] <- as.vector(apply(comparison_diff_draws, c(1, 2), mean))

  overall_comparison <- data.frame(
    Study_Arm = arm_ids,
    One_step = apply(apply(one_step_draws, c(2, 3), sum), 1, mean),
    DE_plus_IE = apply(apply(de_plus_ie_draws, c(2, 3), sum), 1, mean),
    Diff = apply(apply(comparison_diff_draws, c(2, 3), sum), 1, mean)
  )
  names(overall_comparison)[2:4] <- paste0(effect_label, c("_one_step", "_DE_plus_IE", "_diff"))

  list(
    de_plus_ie_draws = de_plus_ie_draws,
    diff_draws = comparison_diff_draws,
    table = comparison_table,
    overall = overall_comparison
  )
}

initialize_arm_result <- function(arm_ids) {
  result <- vector("list", length(arm_ids))
  names(result) <- paste0("Arm_", arm_ids)
  result
}

initialize_effect_draws <- function(year_values, arm_ids, n_draws) {
  array(
    0,
    dim = c(length(year_values), length(arm_ids), n_draws),
    dimnames = list(
      as.character(year_values),
      paste0("Arm_", arm_ids),
      paste0("Draw_", seq_len(n_draws))
    )
  )
}

compute_reference_outcomes <- function(reference_profile,
                                       fit_object,
                                       beta_draws,
                                       treatment_names,
                                       obs_noise_draws = NULL) {
  reference_site_outcomes <- compute_regime_site_outcomes(
    regime_profile = reference_profile,
    fit_object = fit_object,
    beta_draws = beta_draws,
    treatment_names = treatment_names,
    obs_noise_draws = obs_noise_draws
  )

  list(
    site_outcomes = reference_site_outcomes,
    mean = average_site_outcomes(reference_site_outcomes)
  )
}

compute_arm_effects <- function(arm_profiles,
                                arm_ids,
                                reference_mean,
                                fit_object,
                                beta_draws,
                                treatment_names,
                                ie_names,
                                obs_noise_draws = NULL) {
  de_site_outcomes <- initialize_arm_result(arm_ids)
  de_mean_list <- initialize_arm_result(arm_ids)
  arm_site_outcomes <- initialize_arm_result(arm_ids)
  arm_mean_list <- initialize_arm_result(arm_ids)

  de_draws <- initialize_effect_draws(year_values, arm_ids, nrow(beta_draws))
  ie_draws <- initialize_effect_draws(year_values, arm_ids, nrow(beta_draws))
  ate_draws <- initialize_effect_draws(year_values, arm_ids, nrow(beta_draws))

  for (arm_index in seq_along(arm_ids)) {
    arm_id <- arm_ids[arm_index]
    arm_profile <- arm_profiles[arm_profiles$Study_Arm == arm_id, c("Year", treatment_names), drop = FALSE]
    de_profile <- make_de_only_profile(arm_profile, ie_names)

    de_site_outcomes[[arm_index]] <- compute_regime_site_outcomes(
      regime_profile = de_profile,
      fit_object = fit_object,
      beta_draws = beta_draws,
      treatment_names = treatment_names,
      obs_noise_draws = obs_noise_draws
    )
    de_mean_list[[arm_index]] <- average_site_outcomes(de_site_outcomes[[arm_index]])

    arm_site_outcomes[[arm_index]] <- compute_regime_site_outcomes(
      regime_profile = arm_profile,
      fit_object = fit_object,
      beta_draws = beta_draws,
      treatment_names = treatment_names,
      obs_noise_draws = obs_noise_draws
    )
    arm_mean_list[[arm_index]] <- average_site_outcomes(arm_site_outcomes[[arm_index]])

    de_draws[, arm_index, ] <- de_mean_list[[arm_index]] - reference_mean
    ie_draws[, arm_index, ] <- arm_mean_list[[arm_index]] - de_mean_list[[arm_index]]
    ate_draws[, arm_index, ] <- arm_mean_list[[arm_index]] - reference_mean
  }

  list(
    de_site_outcomes = de_site_outcomes,
    de_mean_list = de_mean_list,
    arm_site_outcomes = arm_site_outcomes,
    arm_mean_list = arm_mean_list,
    DE_draws = de_draws,
    IE_draws = ie_draws,
    ATE_draws = ate_draws
  )
}

collect_effect_outputs <- function(effect_draws, arm_ids, year_values) {
  de_results <- summarise_effect_draws(effect_draws$DE_draws, arm_ids, year_values, "DE")
  ie_results <- summarise_effect_draws(effect_draws$IE_draws, arm_ids, year_values, "IE")
  ate_results <- summarise_effect_draws(effect_draws$ATE_draws, arm_ids, year_values, "ATE")
  ate_comparison <- build_decomposition_comparison(
    one_step_draws = effect_draws$ATE_draws,
    de_draws = effect_draws$DE_draws,
    ie_draws = effect_draws$IE_draws,
    arm_ids = arm_ids,
    year_values = year_values,
    effect_label = "ATE"
  )

  c(
    effect_draws,
    list(
      DE_summary = de_results$summary,
      DE_lower = de_results$lower,
      DE_upper = de_results$upper,
      DE_overall = de_results$overall,
      DE_overall_summary = de_results$overall_summary,
      DE_table = de_results$table,
      IE_summary = ie_results$summary,
      IE_lower = ie_results$lower,
      IE_upper = ie_results$upper,
      IE_overall = ie_results$overall,
      IE_overall_summary = ie_results$overall_summary,
      IE_table = ie_results$table,
      ATE_summary = ate_results$summary,
      ATE_lower = ate_results$lower,
      ATE_upper = ate_results$upper,
      ATE_overall = ate_results$overall,
      ATE_overall_summary = ate_results$overall_summary,
      ATE_table = ate_results$table,
      ATE_from_DE_IE_draws = ate_comparison$de_plus_ie_draws,
      ATE_comparison_diff_draws = ate_comparison$diff_draws,
      ATE_comparison_table = ate_comparison$table,
      ATE_comparison_overall = ate_comparison$overall
    )
  )
}

arm_profiles <- extract_arm_profiles(Score_Data, treatment_vars)
year_values <- sort(unique(arm_profiles$Year))
reference_profile <- make_reference_regime(year_values)
beta_draws <- sample_beta_draws(CV_Ranalysis)
obs_noise_draws <- sample_obs_noise_draws(CV_Ranalysis, nrow(beta_draws))
arm_ids <- sort(unique(arm_profiles$Study_Arm))

reference_results <- compute_reference_outcomes(
  reference_profile = reference_profile,
  fit_object = CV_Ranalysis,
  beta_draws = beta_draws,
  treatment_names = treatment_vars
)
reference_site_outcomes <- reference_results$site_outcomes
reference_mean <- reference_results$mean

reference_results_noise <- compute_reference_outcomes(
  reference_profile = reference_profile,
  fit_object = CV_Ranalysis,
  beta_draws = beta_draws,
  treatment_names = treatment_vars,
  obs_noise_draws = obs_noise_draws
)
reference_site_outcomes_noise <- reference_results_noise$site_outcomes
reference_mean_noise <- reference_results_noise$mean

latent_outputs <- collect_effect_outputs(
  effect_draws = compute_arm_effects(
    arm_profiles = arm_profiles,
    arm_ids = arm_ids,
    reference_mean = reference_mean,
    fit_object = CV_Ranalysis,
    beta_draws = beta_draws,
    treatment_names = treatment_vars,
    ie_names = indirect_effect_vars
  ),
  arm_ids = arm_ids,
  year_values = year_values
)

noise_outputs <- collect_effect_outputs(
  effect_draws = compute_arm_effects(
    arm_profiles = arm_profiles,
    arm_ids = arm_ids,
    reference_mean = reference_mean_noise,
    fit_object = CV_Ranalysis,
    beta_draws = beta_draws,
    treatment_names = treatment_vars,
    ie_names = indirect_effect_vars,
    obs_noise_draws = obs_noise_draws
  ),
  arm_ids = arm_ids,
  year_values = year_values
)

invisible(list2env(latent_outputs, environment()))
invisible(list2env(
  stats::setNames(noise_outputs, paste0(names(noise_outputs), "_noise")),
  environment()
))

save(
  arm_profiles,
  reference_profile,
  reference_site_outcomes,
  reference_mean,
  reference_site_outcomes_noise,
  reference_mean_noise,
  de_site_outcomes,
  de_mean_list,
  arm_site_outcomes,
  arm_mean_list,
  de_site_outcomes_noise,
  de_mean_list_noise,
  arm_site_outcomes_noise,
  arm_mean_list_noise,
  obs_noise_draws,
  DE_draws,
  DE_summary,
  DE_lower,
  DE_upper,
  DE_overall,
  DE_overall_summary,
  IE_draws,
  IE_summary,
  IE_lower,
  IE_upper,
  IE_overall,
  IE_overall_summary,
  ATE_draws,
  ATE_summary,
  ATE_lower,
  ATE_upper,
  ATE_overall,
  ATE_overall_summary,
  ATE_from_DE_IE_draws,
  ATE_comparison_diff_draws,
  ATE_comparison_table,
  ATE_comparison_overall,
  DE_draws_noise,
  DE_summary_noise,
  DE_lower_noise,
  DE_upper_noise,
  DE_overall_noise,
  DE_overall_summary_noise,
  IE_draws_noise,
  IE_summary_noise,
  IE_lower_noise,
  IE_upper_noise,
  IE_overall_noise,
  IE_overall_summary_noise,
  ATE_draws_noise,
  ATE_summary_noise,
  ATE_lower_noise,
  ATE_upper_noise,
  ATE_overall_noise,
  ATE_overall_summary_noise,
  ATE_from_DE_IE_draws_noise,
  ATE_comparison_diff_draws_noise,
  ATE_comparison_table_noise,
  ATE_comparison_overall_noise,
  DE_table,
  IE_table,
  ATE_table,
  DE_table_noise,
  IE_table_noise,
  ATE_table_noise,
  file = file.path(output_dir, "JSTVC_ATE_by_arm_vs_zero_regime.RData")
)

print(DE_table)
print(DE_overall_summary)
print(IE_table)
print(IE_overall_summary)
print(ATE_table)
print(ATE_overall_summary)
print(ATE_comparison_table)
print(ATE_comparison_overall)
print(DE_table_noise)
print(DE_overall_summary_noise)
print(IE_table_noise)
print(IE_overall_summary_noise)
print(ATE_table_noise)
print(ATE_overall_summary_noise)
print(ATE_comparison_table_noise)
print(ATE_comparison_overall_noise)
