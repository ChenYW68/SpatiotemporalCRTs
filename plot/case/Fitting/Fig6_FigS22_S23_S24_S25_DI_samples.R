# rm(list = ls())
source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
load("./data/Tanzania_Score_Data_r.RData")
source(normalizePath("./JSTVC/R/util.R"))
standard_treatment_vars <- c(
  "Intercept",
  paste0("CWT_", 1:4),
  paste0("SBT_", 1:4),
  "IEt.CWT",
  "IEt.SBT"
)

inverse_cloglog <- function(eta) {
  1 - exp(-exp(eta))
}

cloglog_transform <- function(x) {
  log(-log(1 - x))
}

compute_spatiotemporal_ate_draws <- function(score_data,
                                             fit_object,
                                             treatment_names = standard_treatment_vars,
                                             transform_fn = inverse_cloglog,
                                             beta_sampler = c("vb", "mcmc"),
                                             n_draws) {
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
      random_effect_draw_indices <- sample_random_effect_draw_indices(fit_object, n_draws = n_draws)
      obs_sigma_sq_draws <- sample_obs_sigma_sq_draws_vb(fit_object, n_draws = n_draws)
      sample_beta_draws_vb(fit_object, n_draws = n_draws)
    },
    mcmc = {
      random_effect_draws_override <- extract_random_effect_draws_mcmc(fit_object, n_draws = n_draws)
      obs_sigma_sq_draws <- sample_obs_sigma_sq_draws_mcmc(fit_object, n_draws = n_draws)
      sample_beta_draws_mcmc(fit_object, n_draws = n_draws)
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
    dimnames = list(as.character(year_values), paste0("Arm_", arm_ids), paste0("Draw_", seq_len(nrow(beta_draws))))
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

  ate_draws
}

compute_fixed_only_ate_draws <- function(score_data,
                                         beta_draws,
                                         treatment_names = standard_treatment_vars,
                                         transform_fn = inverse_cloglog) {
  arm_profiles <- extract_arm_profiles(score_data, treatment_names)
  year_values <- sort(unique(arm_profiles$Year))
  arm_ids <- sort(unique(arm_profiles$Study_Arm))
  reference_profile <- make_reference_regime(year_values, treatment_names)

  ate_draws <- array(
    0,
    dim = c(length(year_values), length(arm_ids), nrow(beta_draws)),
    dimnames = list(as.character(year_values), paste0("Arm_", arm_ids), paste0("Draw_", seq_len(nrow(beta_draws))))
  )

  for (arm_index in seq_along(arm_ids)) {
    arm_profile <- arm_profiles[arm_profiles$Study_Arm == arm_ids[arm_index], c("Year", treatment_names), drop = FALSE]
    for (year_index in seq_along(year_values)) {
      year_value <- year_values[year_index]
      arm_x <- as.numeric(arm_profile[arm_profile$Year == year_value, treatment_names, drop = FALSE][1, ])
      ref_x <- as.numeric(reference_profile[reference_profile$Year == year_value, treatment_names, drop = FALSE][1, ])
      eta_arm <- as.vector(beta_draws %*% arm_x)
      eta_ref <- as.vector(beta_draws %*% ref_x)
      ate_draws[year_index, arm_index, ] <- transform_fn(eta_arm) - transform_fn(eta_ref)
    }
  }

  ate_draws
}

summarise_pairwise_ate_differences <- function(ate_draws) {
  cumulative_ate <- apply(ate_draws, c(2, 3), sum)
  out <- NULL
  n_arms <- nrow(cumulative_ate)

  for (i in seq_len(n_arms - 1)) {
    for (j in seq.int(i + 1, n_arms)) {
      diff_draws <- cumulative_ate[i, ] - cumulative_ate[j, ]
      out <- rbind(out, data.frame(k1 = i, k2 = j, mu = mean(diff_draws), sigma.sq = stats::var(diff_draws)))
    }
  }

  out
}

summarise_pairwise_ate_differences_with_samples <- function(ate_draws) {
  cumulative_ate <- apply(ate_draws, c(2, 3), sum)
  out <- NULL
  n_arms <- nrow(cumulative_ate)
  pair_draws <- list()

  for (i in seq_len(n_arms - 1)) {
    for (j in seq.int(i + 1, n_arms)) {
      diff_draws <- cumulative_ate[i, ] - cumulative_ate[j, ]
      pair_draws[[paste0(i, "_", j)]] <- diff_draws
      out <- rbind(
        out,
        data.frame(
          k1 = i,
          k2 = j,
          mu = mean(diff_draws),
          sigma.sq = stats::var(diff_draws),
          q025 = as.numeric(stats::quantile(diff_draws, 0.025)),
          q975 = as.numeric(stats::quantile(diff_draws, 0.975))
        )
      )
    }
  }

  attr(out, "draws") <- cumulative_ate
  attr(out, "pair_draws") <- pair_draws
  out
}

attach_between_method_quantiles <- function(df_ref, df_cmp) {
  pair_draws_ref <- attr(df_ref, "pair_draws")
  pair_draws_cmp <- attr(df_cmp, "pair_draws")

  if (is.null(pair_draws_ref) || is.null(pair_draws_cmp)) {
    stop("Pairwise summaries must carry paired draws for empirical interval comparison.")
  }

  between_q025 <- numeric(nrow(df_ref))
  between_q975 <- numeric(nrow(df_ref))
  for (row_index in seq_len(nrow(df_ref))) {
    pair_key <- paste0(df_ref$k1[row_index], "_", df_ref$k2[row_index])
    diff_draws <- pair_draws_ref[[pair_key]] - pair_draws_cmp[[pair_key]]
    between_q025[row_index] <- as.numeric(stats::quantile(diff_draws, 0.025))
    between_q975[row_index] <- as.numeric(stats::quantile(diff_draws, 0.975))
  }

  df_ref$between.q025 <- between_q025
  df_ref$between.q975 <- between_q975
  df_cmp$between.q025 <- between_q025
  df_cmp$between.q975 <- between_q975

  list(df_ref = df_ref, df_cmp = df_cmp)
}


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# case = 1 for generating Fig 6 and Fig S22
# case = 2 for generating Fig S24
# case = 3 for generating Fig S25
# case = 4 for generating Fig S23

#--------------------------------------------------------------------
#--------------------------------------------------------------------
N <- 5e3

set.seed(1234)
# Case 1: using combined data
if(case == 1)
{
  Score_Data <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
  fig.name    <- "./figure/Fig6_Dist_ATE_all.pdf"
  figS.name   <- "./figure/FigS22_Dist_ATE_all.pdf"

  #VB
  load(paste0("./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
  vb_fit <- CV_Ranalysis

  #MCMC
  load(paste0("./result/case/JSTVC_Xi_IEt_MCMC_pow_decay_loglog_11.RData"))
  mcmc_fit <- CV_Ranalysis

}

# Case 2: Uing Kenya data only
if(case == 2)
{
  Score_Data <- Kenya_Score_Data
  figS.name   <- "./figure/FigS24_Dist_ATE_Kenya.pdf"

  #VB
  load(paste0("./result/case/Kenya_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
  vb_fit <- CV_Ranalysis
  #MCMC
  load(paste0("./result/case/Kenya_JSTVC_Xi_IEt_MCMC_pow_decay_loglog_11.RData"))
  mcmc_fit <- CV_Ranalysis
}

# Case 3: Uing Tanzania data only
if(case == 3)
{
  Score_Data <- Tanzania_Score_Data
  figS.name   <- "./figure/FigS25_Dist_ATE_Tanzania.pdf"

  #VB
  load(paste0("./result/case/Tanzania_JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
  vb_fit <- CV_Ranalysis

  load(paste0("./result/case/Tanzania_JSTVC_Xi_IEt_MCMC_pow_decay_loglog_11.RData"))
  mcmc_fit <- CV_Ranalysis
}

# case 4
if(case == 4)
{
  Score_Data <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
  figS.name   <- "./figure/FigS23_Dist_ATE_JSTVC_and_JSTVC_xi.pdf"

  #VB
  load(paste0("./result/case/JSTVC_Xi_IEt_VB_pow_decay_loglog_11.RData"))
  vb_fit <- CV_Ranalysis

  # INLA posterior beta draws under cloglog scale
  X_vars <- standard_treatment_vars
  Score_Data_cloglog <- as.data.frame(Score_Data)
  Score_Data_cloglog$Prevalence <- cloglog_transform(Score_Data_cloglog$Prevalence)
  JSTVC_xi <- paste0("Prevalence ~ -1 +  ", paste0(X_vars, collapse = " + "))
  JSTVC_xi.fit  <- inla(as.formula(JSTVC_xi),
                        family = "gaussian",
                        data   = Score_Data_cloglog,
                        verbose = FALSE,
                        control.predictor = list(compute = TRUE),
                        control.compute = list(
                          config = TRUE   # REQUIRED for posterior sampling
                        ))
  post.samp <- inla.posterior.sample(N, JSTVC_xi.fit)
  beta_samples <- inla.posterior.sample.eval(
    function(...) {
      c(
        `Intercept`,
        `CWT_1`,
        `CWT_2`,
        `CWT_3`,
        `CWT_4`,
        `SBT_1`,
        `SBT_2`,
        `SBT_3`,
        `SBT_4`,
        `IEt.CWT`,
        `IEt.SBT`
      )
    },
    post.samp
  )
  beta.sample.mcmc <- t(beta_samples)
}

if(case %in% c(1, 2, 3)) {
  JSTVC.ATE <- compute_spatiotemporal_ate_draws(
    score_data = Score_Data,
    fit_object = vb_fit,
    treatment_names = standard_treatment_vars,
    transform_fn = inverse_cloglog,
    beta_sampler = "vb",
    n_draws = N
  )
  DI <- summarise_pairwise_ate_differences_with_samples(JSTVC.ATE)

  JSTVC.ATE.mcmc <- compute_spatiotemporal_ate_draws(
    score_data = Score_Data,
    fit_object = mcmc_fit,
    treatment_names = standard_treatment_vars,
    transform_fn = inverse_cloglog,
    beta_sampler = "mcmc",
    n_draws = N
  )
  DI.mcmc <- summarise_pairwise_ate_differences_with_samples(JSTVC.ATE.mcmc)
  di_pair <- attach_between_method_quantiles(DI.mcmc, DI)
  DI.mcmc <- di_pair$df_ref
  DI <- di_pair$df_cmp
} else if(case == 4) {
  JSTVC.ATE <- compute_spatiotemporal_ate_draws(
    score_data = Score_Data,
    fit_object = vb_fit,
    treatment_names = standard_treatment_vars,
    transform_fn = inverse_cloglog,
    beta_sampler = "vb",
    n_draws = N
  )
  DI <- summarise_pairwise_ate_differences_with_samples(JSTVC.ATE)

  JSTVC.ATE.mcmc <- compute_fixed_only_ate_draws(
    score_data = Score_Data,
    beta_draws = beta.sample.mcmc,
    treatment_names = standard_treatment_vars,
    transform_fn = inverse_cloglog
  )
  DI.mcmc <- summarise_pairwise_ate_differences_with_samples(JSTVC.ATE.mcmc)
  di_pair <- attach_between_method_quantiles(DI.mcmc, DI)
  DI.mcmc <- di_pair$df_ref
  DI <- di_pair$df_cmp
} else {
  JSTVC.DE <- JSTVC.IE <- JSTVC.DE.mcmc <- JSTVC.IE.mcmc <- array(NA, dim = c(5, 6, N),
                                                                  dimnames = list(paste0("201", 1:5), 1:6, c(1:N)))
  for(k in 1:6){
    Da <- Score_Data[Score_Data$Study_Arm %in% c(k), c(2, 9:16, 19:20)]
    setDF(Da)
    for(i in 1:5){
      X <- Da[Da$Year == (2010 + i), ]
      JSTVC.DE[i, k,] <- colMeans(as.matrix(X[, 2:9]) %*% t(beta.sample[, 1:8]))
      JSTVC.IE[i, k,] <- colMeans(as.matrix(X[, 10:11]) %*% t(beta.sample[, 9:10]))
      JSTVC.DE.mcmc[i, k,] <- colMeans(as.matrix(X[, 2:9]) %*% t(beta.sample.mcmc[, 1:8]))
      JSTVC.IE.mcmc[i, k,] <- colMeans(as.matrix(X[, 10:11]) %*% t(beta.sample.mcmc[, 9:10]))
    }
  }
  JSTVC.ATE <- JSTVC.DE + JSTVC.IE
  ATE <- apply(JSTVC.ATE, 2:3, sum)
  DI <- NULL
  DI.draws <- list()
  for(i in 1:5)
  {
    for(j in (i + 1):6){
      diff_draws <- ATE[i,] - ATE[j,]
      DI.draws[[paste0(i, "_", j)]] <- diff_draws
      DI <- rbind(DI, data.frame(k1 = i, k2 = j,
                                 mu = mean(diff_draws),
                                 sigma.sq = var(diff_draws),
                                 q025 = as.numeric(stats::quantile(diff_draws, 0.025)),
                                 q975 = as.numeric(stats::quantile(diff_draws, 0.975))))
    }
  }
  attr(DI, "pair_draws") <- DI.draws

  JSTVC.ATE.mcmc <- JSTVC.DE.mcmc + JSTVC.IE.mcmc
  ATE.mcmc <- apply(JSTVC.ATE.mcmc, 2:3, sum)
  DI.mcmc <- NULL
  DI.mcmc.draws <- list()
  for(i in 1:5)
  {
    for(j in (i + 1):6){
      diff_draws <- ATE.mcmc[i,] - ATE.mcmc[j,]
      DI.mcmc.draws[[paste0(i, "_", j)]] <- diff_draws
      DI.mcmc <- rbind(DI.mcmc, data.frame(k1 = i, k2 = j,
                                           mu = mean(diff_draws),
                                           sigma.sq = var(diff_draws),
                                           q025 = as.numeric(stats::quantile(diff_draws, 0.025)),
                                           q975 = as.numeric(stats::quantile(diff_draws, 0.975))))
    }
  }
  attr(DI.mcmc, "pair_draws") <- DI.mcmc.draws
  di_pair <- attach_between_method_quantiles(DI.mcmc, DI)
  DI.mcmc <- di_pair$df_ref
  DI <- di_pair$df_cmp
}


lab <- c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)",
         "(G)", "(H)", "(J)", "(K)", "(L)", "(M)",
         "(N)", "(O)", "(P)")

Ind <- c(1, 14, 15)#c(2, 13, 14) #c(3, 14, 15)
if(case == 1){
  df      <- DI[Ind, ]
  df.mcmc <- DI.mcmc[Ind, ]

  # Open a PDF device to save the plots (each page will have one figure)
  pdf(file  = fig.name,
      width = 21,
      height = 6)

  par(mar = c(3, 3.5, 2, 1), mgp = c(2, 0.8, 0))
  par(mfrow = c(1, 3),
      cex = 1.2,
      cex.axis = 1.4,
      cex.lab = 1.5,
      cex.main = 1.4,
      lwd = 0.5)

  ylab  <- "Density"
  # Loop over each row to produce a plot
  for (i in 1:nrow(df)) {
    if(i != 1){
      ylab  <- ""
    }
    plot_two_normal_distributions_samples(mean1 = df.mcmc$mu[i],
                                          sd1   = sqrt(df.mcmc$sigma.sq[i]),
                                          mean2 = df$mu[i],
                                          sd2   = sqrt(df$sigma.sq[i]),
                                          ci1_text = c(df.mcmc$q025[i], df.mcmc$q975[i]),
                                          ci2_text = c(df$q025[i], df$q975[i]),
                                          ci_between_text = c(df.mcmc$between.q025[i], df.mcmc$between.q975[i]),
                                          main  = paste(lab[i], "Arm ", df$k1[i], "vs", df$k2[i]),
                                          hjust = 0.2,
                                          ylab  = ylab,
                                          bar.cex  = 1,
                                          text.cex = 1,
                                          x.adjust = 0.001,
                                          y.pos = c(0.895, 0.84, 0.72, 0.60))

    # if(i == 1){
    #   mtext(paste(lab[i], "Arm ", df$k1[i], "vs", df$k2[i]), side = 3, line = 0.51, adj = 0, cex = 1.8)
    # }else{
    #   mtext(paste(lab[i], "Arm ", df$k1[i], "vs", df$k2[i]), side = 3, line = 0.51, adj = 0, cex = 1.8)
    # }

  }

  # Close the PDF device to save the file
  dev.off()
  # if(i == 1){
  cat(sprintf(
    "Arm %s vs %s: DI_12 = %.3f [%.3f, %.3f]\n",
    df$k1[1],
    df$k2[1],
    df$mu[1],
    df$q025[1],
    df$q975[1]
  ))
  # }
}



if(case == 1){
  pdf(file  = figS.name,
      width = 18,
      height = 20)
  par(mar = c(3, 3.5, 2.5, 1), mgp = c(2, 0.8, 0))
  par(mfrow = c(4, 3),
      cex = 1.2,
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.3,
      lwd = 0.5)
  df      <- DI[-Ind, ]
  df.mcmc <- DI.mcmc[-Ind, ]
  for (i in 1:nrow(df)) {
    ylab <- "Density"
    if(i %nin% c(1, 4, 7, 10, 13)){
      ylab  <- ""
    }
    plot_two_normal_distributions_samples(mean1 = df.mcmc$mu[i],
                                          sd1   = sqrt(df.mcmc$sigma.sq[i]),
                                          mean2 = df$mu[i],
                                          sd2   = sqrt(df$sigma.sq[i]),
                                          ci1_text = c(df.mcmc$q025[i], df.mcmc$q975[i]),
                                          ci2_text = c(df$q025[i], df$q975[i]),
                                          ci_between_text = c(df.mcmc$between.q025[i], df.mcmc$between.q975[i]),
                                          main  = paste(lab[i], "Arm ", df$k1[i], "vs", df$k2[i]),
                                          hjust = 0.2,
                                          ylab  = ylab,
                                          bar.cex  = 0.8,
                                          text.cex = 0.75,
                                          x.adjust = 0.001,
                                          y.pos = c(0.895, 0.84, 0.70, 0.56))
  }
  dev.off()
}
if(case%in% c(2, 3)){
  pdf(file  = figS.name,
      width = 20,
      height = 28)
  par(mar = c(3.5, 4, 3, 2), mgp = c(2, 0.8, 0))
  par(mfrow = c(5, 3),
      cex = 1.3,
      cex.axis = 1.3,
      cex.lab = 1.5,
      cex.main = 1.2,
      lwd = 1)
  df      <- DI
  df.mcmc <- DI.mcmc
  for (i in 1:nrow(df)) {
    ylab <- "Density"
    if(i %nin% c(1, 4, 7, 10, 13)){
      ylab  <- ""
    }
    plot_two_normal_distributions_samples(mean1 = df.mcmc$mu[i],
                                          sd1   = sqrt(df.mcmc$sigma.sq[i]),
                                          mean2 = df$mu[i],
                                          sd2   = sqrt(df$sigma.sq[i]),
                                          ci1_text = c(df.mcmc$q025[i], df.mcmc$q975[i]),
                                          ci2_text = c(df$q025[i], df$q975[i]),
                                          ci_between_text = c(df.mcmc$between.q025[i], df.mcmc$between.q975[i]),
                                          main  = paste(lab[i], "Arm ", df$k1[i], "vs", df$k2[i]),
                                          hjust = 0.2,
                                          ylab  = ylab,
                                          bar.cex  = 0.8,
                                          text.cex = 0.75,
                                          x.adjust = 0.001,
                                          y.pos = c(0.895, 0.84, 0.70, 0.56))
  }
  dev.off()
}

if(case == 4){
  pdf(file  = figS.name,
      width = 20,
      height = 28)
  par(mar = c(3.5, 4, 3, 2), mgp = c(2, 0.8, 0))
  par(mfrow = c(5, 3),
      cex = 1.3,
      cex.axis = 1.3,
      cex.lab = 1.5,
      cex.main = 1.2,
      lwd = 1)
  df      <- DI
  df.mcmc <- DI.mcmc
  for (i in 1:nrow(df)) {
    ylab <- "Density"
    if(i %nin% c(1, 4, 7, 10, 13)){
      ylab  <- ""
    }
    plot_two_normal_distributions_samples(mean1 = df.mcmc$mu[i],
                                          sd1   = sqrt(df.mcmc$sigma.sq[i]),
                                          mean2 = df$mu[i],
                                          sd2   = sqrt(df$sigma.sq[i]),
                                          ci1_text = c(df.mcmc$q025[i], df.mcmc$q975[i]),
                                          ci2_text = c(df$q025[i], df$q975[i]),
                                          ci_between_text = c(df.mcmc$between.q025[i], df.mcmc$between.q975[i]),
                                          main  = paste(lab[i], "Arm ", df$k1[i], "vs", df$k2[i]),
                                          labels   = c(TeX("JSTVC$_{-\\xi}$"), "JSTVC"),
                                          hjust = 0,
                                          ylab  = ylab,
                                          bar.cex  = 0.8,
                                          text.cex = 0.75,
                                          x.adjust = 0.001,
                                          y.pos = c(0.890, 0.83, 0.70, 0.55),
                                          show.CI.label = FALSE)
  }
  dev.off()
}
