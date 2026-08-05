#' Calculate performance summaries
#'
#' @param estimates Matrix or data frame of estimates.
#' @param true_values Vector of true values.
calc_performance <- function(estimates, true_values) {
  # estimates: matrix or data frame of size (n_rounds x n_coef)
  # true_values: vector of true coefficient values (length = n_coef)

  # Check dimensions
  if (ncol(estimates) != length(true_values)) {
    stop("Number of columns in estimates must equal length of true_values.")
  }

  # Compute summary statistics
  bias <- colMeans(estimates) - true_values
  sd <- apply(estimates, 2, sd)
  rmse <- (colMeans((estimates - matrix(true_values,
                                        nrow = nrow(estimates),
                                        ncol = length(true_values),
                                        byrow = TRUE))^2))

  # Combine into a data frame
  Details <- data.frame(
    Coefficient = paste0("Beta", seq_along(true_values)),
    Bias = round(bias, 5),
    SD = round(sd, 5),
    MSE = round(rmse, 5)
  )

  return(list(Details = Details))
}

inverse_cloglog <- function(eta) {
  1 - exp(-exp(eta))
}

cloglog_transform <- function(x) {
  log(-log(1 - x))
}

extract_arm_profiles <- function(score_data, treatment_names) {
  score_data |>
    dplyr::select(Study_Arm, Year, dplyr::all_of(treatment_names)) |>
    dplyr::group_by(Study_Arm, Year) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(treatment_names), ~ .x[1]), .groups = "drop") |>
    dplyr::arrange(Study_Arm, Year)
}

make_reference_regime <- function(years, treatment_names) {
  reference_profile <- as.data.frame(matrix(0, nrow = length(years), ncol = length(treatment_names)))
  names(reference_profile) <- treatment_names
  reference_profile$Year <- years
  reference_profile$Intercept <- 1
  reference_profile[, c("Year", treatment_names), drop = FALSE]
}

sample_beta_draws_vb <- function(fit_object, n_draws) {
  sigma <- fit_object$update.Para.List[[1]]$beta$sigma.sq
  if (length(dim(sigma)) != 2) {
    sigma <- as.matrix(sigma)
  }

  mvnfast::rmvn(
    n = n_draws,
    mu = as.vector(fit_object$update.Para.List[[1]]$beta$mu.beta),
    sigma = sigma
  )
}

sample_beta_draws_mcmc <- function(fit_object, n_draws) {
  beta_chain <- fit_object$Process.monitoring$VB.EnKS.Iter |>
    dplyr::select(dplyr::starts_with("pub."))

  if (nrow(beta_chain) < n_draws) {
    stop("The MCMC chain does not contain enough posterior draws.")
  }

  as.matrix(utils::tail(beta_chain, n_draws))
}

sample_obs_sigma_sq_draws_vb <- function(fit_object, n_draws) {
  region_names <- names(fit_object$data)
  stats::setNames(
    lapply(seq_along(region_names), function(region_index) {
      sigma_info <- fit_object$update.Para.List[[region_index]]$obs.sigma.sq
      if (!is.null(sigma_info$a) && !is.null(sigma_info$b)) {
        return(1 / stats::rgamma(n_draws, shape = sigma_info$a, rate = sigma_info$b))
      }

      rep(sigma_info$mu.sigma.sq, n_draws)
    }),
    region_names
  )
}

sample_obs_sigma_sq_draws_mcmc <- function(fit_object, n_draws) {
  sigma_chain <- fit_object$Process.monitoring$VB.EnKS.Iter
  region_names <- names(fit_object$data)
  sigma_draws <- vector("list", length(region_names))
  names(sigma_draws) <- region_names

  for (region_index in seq_along(region_names)) {
    sigma_col <- paste0(region_names[region_index], ".obs.sigma.sq")
    if (!sigma_col %in% names(sigma_chain)) {
      stop(paste0("The MCMC chain does not contain observation variance draws for ", region_names[region_index], "."))
    }

    sigma_values <- utils::tail(sigma_chain[[sigma_col]], n_draws)
    if (length(sigma_values) < n_draws) {
      stop("The MCMC chain does not contain enough observation variance draws.")
    }
    sigma_draws[[region_index]] <- sigma_values
  }

  sigma_draws
}

sample_random_effect_draw_indices <- function(fit_object, n_draws) {
  region_names <- names(fit_object$data)
  if (length(region_names) == 0) {
    stop("The fit object does not contain any regions for random-effect sampling.")
  }

  random_effect_array <- fit_object$fitted.slope.Hv.Zg.ens[[paste0(region_names[1], ".Intercept")]]
  available_draws <- dim(random_effect_array)[3]

  if (is.null(available_draws) || available_draws == 0) {
    stop("The fit object does not contain spatiotemporal random-effect draws.")
  }

  sample.int(available_draws, size = n_draws, replace = available_draws < n_draws)
}

extract_random_effect_draws_mcmc <- function(fit_object, n_draws) {
  random_effect_history <- fit_object$Process.monitoring$MCMC.random.effect.history
  region_names <- names(fit_object$data)

  if (is.null(random_effect_history)) {
    stop("The fit object does not contain MCMC random-effect history.")
  }

  random_effect_draws <- vector("list", length(region_names))
  names(random_effect_draws) <- region_names

  for (region_index in seq_along(region_names)) {
    region_name <- paste0(region_names[region_index], ".Intercept")
    if (!region_name %in% names(random_effect_history)) {
      stop(paste0("The MCMC random-effect history does not contain draws for ", region_name, "."))
    }

    region_history <- random_effect_history[[region_name]][
      ,
      fit_object$data[[region_index]]$index.y,
      ,
      drop = FALSE
    ]

    available_draws <- dim(region_history)[3]
    if (is.null(available_draws) || available_draws == 0) {
      stop("The MCMC random-effect history does not contain posterior draws.")
    }
    if (available_draws < n_draws) {
      stop("The MCMC random-effect history does not contain enough posterior draws.")
    }

    random_effect_draws[[region_names[region_index]]] <- region_history[
      ,
      ,
      seq.int(available_draws - n_draws + 1, available_draws),
      drop = FALSE
    ]
  }

  random_effect_draws
}

compute_regime_site_outcomes <- function(regime_profile,
                                         fit_object,
                                         beta_draws,
                                         obs_sigma_sq_draws = NULL,
                                         treatment_names,
                                         transform_fn,
                                         random_effect_draw_indices = NULL,
                                         random_effect_draws_override = NULL) {
  region_names <- names(fit_object$data)
  n_draws <- nrow(beta_draws)
  year_values <- sort(unique(regime_profile$Year))
  n_years <- length(year_values)

  site_outcomes <- vector("list", length(region_names))
  names(site_outcomes) <- region_names

  for (region_index in seq_along(region_names)) {
    region_name <- region_names[region_index]
    random_effect_draws <- random_effect_draws_override[[region_name]]
    if (is.null(random_effect_draws)) {
      random_effect_draws <- fit_object$fitted.slope.Hv.Zg.ens[[paste0(region_name, ".Intercept")]][
        ,
        fit_object$data[[region_index]]$index.y,
        ,
        drop = FALSE
      ]
      if (!is.null(random_effect_draw_indices)) {
        random_effect_draws <- random_effect_draws[, , random_effect_draw_indices, drop = FALSE]
      }
    }
    if (dim(random_effect_draws)[3] != n_draws) {
      stop("Random-effect draws must match the number of posterior draws.")
    }

    sigma_sq_draws <- NULL
    if (!is.null(obs_sigma_sq_draws)) {
      sigma_sq_draws <- obs_sigma_sq_draws[[region_name]]
      if (length(sigma_sq_draws) != n_draws) {
        stop("Observation variance draws must match the number of posterior draws.")
      }
    }

    n_sites <- dim(random_effect_draws)[2]
    site_ids <- colnames(fit_object$data[[region_index]]$Y_ts)
    if (is.null(site_ids)) {
      site_ids <- paste0(region_name, "_Site_", seq_len(n_sites))
    }

    site_outcomes[[region_name]] <- array(
      NA_real_,
      dim = c(n_years, n_sites, n_draws),
      dimnames = list(as.character(year_values), site_ids, paste0("Draw_", seq_len(n_draws)))
    )

    for (year_index in seq_len(n_years)) {
      year_value <- year_values[year_index]
      profile_row <- regime_profile[regime_profile$Year == year_value, treatment_names, drop = FALSE]

      if (nrow(profile_row) != 1) {
        stop("Each regime must provide exactly one treatment profile per year.")
      }

      x_vector <- as.numeric(profile_row[1, ])
      x_beta_draws <- as.vector(beta_draws %*% x_vector)
      random_effect_matrix <- matrix(
        random_effect_draws[year_index, , ],
        nrow = n_sites,
        ncol = n_draws
      )
      eta_draws <- sweep(random_effect_matrix, 2, x_beta_draws, "+")

      if (!is.null(sigma_sq_draws)) {
        noise_matrix <- matrix(
          stats::rnorm(n_sites * n_draws),
          nrow = n_sites,
          ncol = n_draws
        )
        noise_matrix <- sweep(noise_matrix, 2, sqrt(sigma_sq_draws), "*")
        eta_draws <- eta_draws + noise_matrix
      }

      site_outcomes[[region_name]][year_index, , ] <- transform_fn(eta_draws)
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
#' Compute empirical coverage indicators
#'
#' @param beta_mean Matrix of posterior means.
#' @param beta_var Matrix of posterior variances.
#' @param beta_true Vector of true coefficient values.
#' @param conf_level Confidence level used for interval construction.
compute_coverage <- function(beta_mean, beta_var, beta_true, conf_level = 0.95) {
  # beta_mean: R x p matrix of posterior means
  # beta_var: R x p matrix of posterior variances
  # beta_true: vector of true beta values (length p)
  # conf_level: confidence level (default 0.95)

  R <- nrow(beta_mean)
  p <- ncol(beta_mean)

  z <- qnorm(1 - (1 - conf_level) / 2)  # two-sided z

  coverage_matrix <- matrix(FALSE, R, p)

  for (r in 1:R) {
    for (j in 1:p) {
      sd_j <- sqrt(beta_var[r, j])
      lower <- beta_mean[r, j] - z * sd_j
      upper <- beta_mean[r, j] + z * sd_j
      coverage_matrix[r, j] <- (beta_true[j] >= lower && beta_true[j] <= upper)
    }
  }

  return(coverage_matrix)
}

#' Build spatial taper matrices
#'
#' @param cs Spatial tapering range or ranges.
#' @param data Input data object.
spatial.taper <- function(cs, data){
  spTaper <- list()
  if(length(cs) < length(data)){
    cs <- rep(cs, length(data))
  }
  for(py in 1:length(data)){
    spTaper[[py]] <- list()
    for(g in 1:data[[py]]$Grid.infor$summary$res) {
      spTaper[[py]]$tuning[g] <- cs[py]
      spTaper[[py]]$taper[[g]] <- spam::as.spam(fields::Wendland(data[[py]]$Grid.infor$level[[g]]$BAUs.Dist
                                                , theta = spTaper[[py]]$tuning[g]
                                                , dimension = 1, k = 1))
    }
  }

  return(spTaper)
}
# tapering function
#' Evaluate the Wendland taper function
#'
#' @param distance Distance values.
#' @param radius Taper radius.
taper <- function(distance, radius) Wendland(distance, radius, 1, 1)


#' Plot a normal density
#'
#' @param mean Mean of the normal distribution.
#' @param sd Standard deviation of the normal distribution.
#' @param from Lower bound for the x-axis.
#' @param to Upper bound for the x-axis.
#' @param n Number of grid points.
#' @param main Plot title.
#' @param col Line color.
#' @param lwd Line width.
#' @param ... Additional arguments passed to `plot()`.
plot_normal <- function(mean = 0, sd = 1, from = NULL, to = NULL, n = 200,
                        main = "Normal Density", col = "blue", lwd = 2, ...) {
  # Default x-range: ±4 standard deviations around the mean
  if (is.null(from)) from <- mean - 4 * sd
  if (is.null(to))   to <- mean + 4 * sd

  # Generate x values
  x <- seq(from, to, length.out = n)
  # Compute density values
  y <- dnorm(x, mean = mean, sd = sd)

  # Plot
  plot(x, y, type = "l", main = main,
       xlab = "x", ylab = "Density", col = col, lwd = lwd, ...)

  # Add vertical line at the mean
  abline(v = mean, col = "red", lty = 2)
}

#' Plot one normal distribution with interval shading
#'
#' @param mean1 Mean of the distribution.
#' @param sd1 Standard deviation of the distribution.
#' @param xm X-axis limits.
#' @param ym Y-axis upper limit.
#' @param lab.p Panel label.
plot_two_normal_distribution <- function(mean1, sd1, xm, ym, lab.p) {
  # Define the range for x values
  x <- seq(min(mean1 - 4*sd1),
           max(mean1 + 4*sd1),
           length.out = 100)

  # Calculate the density values for both distributions
  y1 <- dnorm(x, mean = mean1, sd = sd1)


  # Create the plot for the first normal distribution
  plot(x, y1, type = 'l', lwd = 2, col = 'gray',
       main = "Baseline",#latex2exp::TeX("$\\beta_0$")
       xlab = "Intercept", ylab = "Density",
       xlim = c(xm[1], xm[2]),
       ylim = c(0, ym))
  abline(v = mean1, col = scales::alpha("gray", alpha = 0.3))
  mtext(lab.p, side = 3, line = 0.51, adj = 0, cex = 2)
  ci1 <- qnorm(c(0.025, 0.975), mean = mean1, sd = sd1)
  polygon(c(ci1[1], seq(ci1[1], ci1[2], length.out = 100), ci1[2]),
          c(0, dnorm(seq(ci1[1], ci1[2], length.out = 100), mean = mean1, sd = sd1), 0),
          col = rgb(0.5, 0.5, 0.5, alpha = 0.2), border = NA)
  # Add the second normal distribution to the plot

  # Add a legend
  legend("topleft", legend = c(paste("N(", Round(mean1, 2), ",", Round(sd1, 2), ")", sep = "")),
         col = c("gray", "red"), lwd = 2)

  # Add a grid for better visibility
  grid()
}

#' Plot two normal distributions for comparison
#'
#' @param mean1 Mean of the first distribution.
#' @param sd1 Standard deviation of the first distribution.
#' @param mean2 Mean of the second distribution.
#' @param sd2 Standard deviation of the second distribution.
#' @param main Plot title.
#' @param xm X-axis limits.
#' @param ym Y-axis upper limit.
#' @param hjust Horizontal justification setting.
#' @param vjust Vertical justification setting.
#' @param lab.p Panel label.
#' @param ylab Y-axis label.
#' @param labels Labels for the two distributions.
#' @param bar.cex Character expansion for the legend.
#' @param text.cex Character expansion for annotations.
#' @param x.adjust Horizontal position adjustment for text.
#' @param y.pos Relative y positions for annotations.
#' @param show.CI.label Logical; whether to show CI endpoint labels.
plot_two_normal_distributions <- function(mean1, sd1,
                                          mean2, sd2,
                                          main,
                                          xm = NULL,
                                          ym = NULL,
                                          hjust = 0,
                                          vjust = 0,
                                          lab.p = NULL,
                                          ylab  = "Density",
                                          labels   = c("MCMC", "VB"),
                                          bar.cex  = 1,
                                          text.cex = 1,
                                          x.adjust = 0.05,
                                          y.pos    = c(0.9, 0.84, 0.72, 0.60),
                                          show.CI.label = TRUE) {
  # Define the range for x values
  x <- seq(min(mean1 - 4*sd1, mean2 - 4*sd2),
           max(mean1 + 4*sd1, mean2 + 4*sd2),
           length.out = 100)

  y1 <- dnorm(x, mean = mean1, sd = sd1)
  y2 <- dnorm(x, mean = mean2, sd = sd2)

  if(is.null(xm)){
    xm <- c(min(x), max(x))
  }
  if(is.null(ym)){
    ym <- max(y1, y2)
  }

  # Create the plot for the first normal distribution
  plot(x, y1,
       type = 'l',
       lwd = 3,
       lty = 5,
       col = 'blue',
       xlab = "Difference in ATEs",
       ylab = ylab,
       xlim = c(xm[1], xm[2]),
       ylim = c(0, ym)
  )
  # abline(v = mean1, col = scales::alpha("blue", alpha = 0.3))
  ci1 <- qnorm(c(0.025, 0.975), mean = mean1, sd = sd1)
  if(!is.null(lab.p)){
    mtext(lab.p, side = 3, line = 0.51, adj = 0, cex = 2)
  }
  lower1 <- mean1 - 1.96 * sd1
  upper1 <- mean1 + 1.96 * sd1

  y_lower1 <- dnorm(lower1, mean = mean1, sd = sd1)
  y_upper1 <- dnorm(upper1, mean = mean1, sd = sd1)

  segments(x0 = lower1, y0 = 0, x1 = lower1, y1 = y_lower1,
           col = scales::alpha('blue', alpha = 0.5), lty = 5, lwd = 3)
  segments(x0 = upper1, y0 = 0, x1 = upper1, y1 = y_upper1,
           col = scales::alpha('blue', alpha = 0.4), lty = 5, lwd = 3)

  polygon(c(ci1[1], seq(ci1[1], ci1[2], length.out = 100), ci1[2]),
          c(0, dnorm(seq(ci1[1], ci1[2], length.out = 100), mean = mean1, sd = sd1), 0),
          col = rgb(0, 0, 1, alpha = 0.2), border = NA)

  # Add the second normal distribution to the plot
  lines(x, y2, lwd = 2, col = 'red')
  # abline(v = mean2, col = scales::alpha("red", alpha = 0.3))


  ci2 <- qnorm(c(0.025, 0.975), mean = mean2, sd = sd2)
  polygon(c(ci2[1], seq(ci2[1], ci2[2], length.out = 100), ci2[2]),
          c(0, dnorm(seq(ci2[1], ci2[2], length.out = 100), mean = mean2, sd = sd2), 0),
          col = rgb(1, 0, 0, alpha = 0.2), border = NA)

  lower2 <- mean2 - 1.96 * sd2
  upper2 <- mean2 + 1.96 * sd2

  y_lower2 <- dnorm(lower2, mean = mean2, sd = sd2)
  y_upper2 <- dnorm(upper2, mean = mean2, sd = sd2)

  segments(x0 = lower2, y0 = 0, x1 = lower2, y1 = y_lower2,
           col = scales::alpha('red', alpha = 0.5), lty = 1, lwd = 3)
  segments(x0 = upper2, y0 = 0, x1 = upper2, y1 = y_upper2,
           col = scales::alpha('red', alpha = 0.4), lty = 1, lwd = 3)

 if(show.CI.label){
   text(max(lower1, lower2), min(y_lower1, y_lower2), labels = sprintf("%.2f", lower1),
        pos = 4, col = "blue")
   text(min(upper1, upper2), min(y_upper1, y_upper2), labels = sprintf("%.2f", upper1),
        pos = 2, col = "blue")


   text(max(lower1, lower2), 1e-1, labels = sprintf("%.2f", lower2),
        pos = 4, col = "red")
   text(min(upper1, upper2), 1e-1, labels = sprintf("%.2f", upper2),
        pos = 2, col = "red")
 }

  posterior_diff <- rnorm(1E4, mean1 - mean2, sqrt(sd1^2 + sd2^2))

  # Calculate the 95% credible interval of the mean difference
  credible_interval <- quantile(posterior_diff, c(0.025, 0.975))


  usr <- par("usr")  # c(xmin, xmax, ymin, ymax)
  legend_x <- usr[1] - 0.02*(usr[2]-usr[1])
  legend_y <- usr[4] + 0.02*(usr[4]-usr[3])

  leg <- legend(#"topleft",
                x = legend_x,
                y = legend_y,
         legend = labels,
         col = c("blue", "red"),
         lwd = c(3, 3),
         lty = c(5, 1),
         ncol = 2,
         bty    = "n",
         cex   = bar.cex,
         seg.len = 3,
         inset = c(-0.01, 0.01))

  y_top <- leg$rect[4]$top


  y_rel <- y.pos

  y_positions <- usr[3] + y_rel * (y_top - usr[3])


  text(
    x = usr[1] + x.adjust,
    y = y_positions[1],
    labels = "95% CIs of ATE difference",
    font = 2,
    col = "black",
    cex = text.cex,
    adj = c(0, 1)
  )


  text(
    x = usr[1] + x.adjust,
    y = y_positions[2],
    labels = as.expression(bquote(
      atop(
        "Between arms using " * .(labels[[1]]) * ":",
        paste("[" * .(Round(lower1, 2)) * ", " * .(Round(upper1, 2)) * "]      ")
      )
    )),#paste0("Between arms using", labels[1], ":\n[", Round(lower1,2), ", ", Round(upper1,2), "]"),
    col = "black",
    cex = text.cex,
    adj = c(0,1)
  )

  text(
    x = usr[1] + x.adjust,
    y = y_positions[3],
    labels =  as.expression(bquote(
      atop(
        "Between arms using " * .(labels[[2]]) * ":",
        "[" * .(Round(lower2, 2)) * ", " * .(Round(upper2, 2)) * "]"
      )
    )),#paste0("Between arms using VB:\n[", Round(lower2,2), ", ", Round(upper2,2), "]"),
    col = "black",
    cex = text.cex,
    adj = c(0,1)
  )


  text(
    x = usr[1] + x.adjust,
    y = y_positions[4],
    labels = as.expression(bquote(atop(
      "Between " * .(labels[[1]]) * " and " * .(labels[[2]]) * ":",
      "[" * .(Round(credible_interval[1], 2)) * ", " * .(Round(credible_interval[2], 2)) * "]"
    ))),
      #paste0("Between", labels[1], "and" , labels[2], ":\n[", Round(credible_interval[1],2), ", ", Round(credible_interval[2],2), "]"),
    col = "black",
    cex = text.cex,
    adj = c(0,1)
  )

  mtext(main, side = 3, line = 0.51, adj = 0, cex = 1.8)
  # Add a grid for better visibility
  # grid()
}

plot_two_normal_distributions_samples <- function(mean1, sd1,
                                                  mean2, sd2,
                                                  ci1_text,
                                                  ci2_text,
                                                  ci_between_text,
                                                  main,
                                                  xm = NULL,
                                                  ym = NULL,
                                                  hjust = 0,
                                                  vjust = 0,
                                                  lab.p = NULL,
                                                  ylab  = "Density",
                                                  labels   = c("MCMC", "VB"),
                                                  bar.cex  = 1,
                                                  text.cex = 1,
                                                  x.adjust = 0.05,
                                                  y.pos    = c(0.9, 0.84, 0.72, 0.60),
                                                  show.CI.label = TRUE) {
  x <- seq(min(mean1 - 4 * sd1, mean2 - 4 * sd2),
           max(mean1 + 4 * sd1, mean2 + 4 * sd2),
           length.out = 100)

  y1 <- dnorm(x, mean = mean1, sd = sd1)
  y2 <- dnorm(x, mean = mean2, sd = sd2)

  if (is.null(xm)) {
    xm <- c(min(x), max(x))
  }
  if (is.null(ym)) {
    ym <- max(y1, y2)
  }

  plot(x, y1,
       type = "l",
       lwd = 3,
       lty = 5,
       col = "blue",
       xlab = "Difference in ATEs",
       ylab = ylab,
       xlim = c(xm[1], xm[2]),
       ylim = c(0, ym))

  ci1_normal <- qnorm(c(0.025, 0.975), mean = mean1, sd = sd1)
  if (!is.null(lab.p)) {
    mtext(lab.p, side = 3, line = 0.51, adj = 0, cex = 2)
  }

  lower1 <- ci1_normal[1]
  upper1 <- ci1_normal[2]
  y_lower1 <- dnorm(lower1, mean = mean1, sd = sd1)
  y_upper1 <- dnorm(upper1, mean = mean1, sd = sd1)

  segments(x0 = lower1, y0 = 0, x1 = lower1, y1 = y_lower1,
           col = scales::alpha("blue", alpha = 0.5), lty = 5, lwd = 3)
  segments(x0 = upper1, y0 = 0, x1 = upper1, y1 = y_upper1,
           col = scales::alpha("blue", alpha = 0.4), lty = 5, lwd = 3)
  polygon(c(ci1_normal[1], seq(ci1_normal[1], ci1_normal[2], length.out = 100), ci1_normal[2]),
          c(0, dnorm(seq(ci1_normal[1], ci1_normal[2], length.out = 100), mean = mean1, sd = sd1), 0),
          col = rgb(0, 0, 1, alpha = 0.2), border = NA)

  lines(x, y2, lwd = 2, col = "red")

  ci2_normal <- qnorm(c(0.025, 0.975), mean = mean2, sd = sd2)
  polygon(c(ci2_normal[1], seq(ci2_normal[1], ci2_normal[2], length.out = 100), ci2_normal[2]),
          c(0, dnorm(seq(ci2_normal[1], ci2_normal[2], length.out = 100), mean = mean2, sd = sd2), 0),
          col = rgb(1, 0, 0, alpha = 0.2), border = NA)

  lower2 <- ci2_normal[1]
  upper2 <- ci2_normal[2]
  y_lower2 <- dnorm(lower2, mean = mean2, sd = sd2)
  y_upper2 <- dnorm(upper2, mean = mean2, sd = sd2)

  segments(x0 = lower2, y0 = 0, x1 = lower2, y1 = y_lower2,
           col = scales::alpha("red", alpha = 0.5), lty = 1, lwd = 3)
  segments(x0 = upper2, y0 = 0, x1 = upper2, y1 = y_upper2,
           col = scales::alpha("red", alpha = 0.4), lty = 1, lwd = 3)

  if (show.CI.label) {
    text(max(lower1, lower2), min(y_lower1, y_lower2), labels = sprintf("%.2f", lower1),
         pos = 4, col = "blue")
    text(min(upper1, upper2), min(y_upper1, y_upper2), labels = sprintf("%.2f", upper1),
         pos = 2, col = "blue")
    text(max(lower1, lower2), 1e-1, labels = sprintf("%.2f", lower2),
         pos = 4, col = "red")
    text(min(upper1, upper2), 1e-1, labels = sprintf("%.2f", upper2),
         pos = 2, col = "red")
  }

  usr <- par("usr")
  legend_x <- usr[1] - 0.02 * (usr[2] - usr[1])
  legend_y <- usr[4] + 0.02 * (usr[4] - usr[3])

  leg <- legend(
    x = legend_x,
    y = legend_y,
    legend = labels,
    col = c("blue", "red"),
    lwd = c(3, 3),
    lty = c(5, 1),
    ncol = 2,
    bty = "n",
    cex = bar.cex,
    seg.len = 3,
    inset = c(-0.01, 0.01)
  )

  y_top <- leg$rect[4]$top
  y_positions <- usr[3] + y.pos * (y_top - usr[3])

  text(x = usr[1] + x.adjust, y = y_positions[1],
       labels = "95% CIs of mean difference",
       font = 2, col = "black", cex = text.cex, adj = c(0, 1))

  text(
    x = usr[1] + x.adjust,
    y = y_positions[2],
    labels = as.expression(bquote(
      atop(
        "Between arms using " * .(labels[[1]]) * ":",
        "[" * .(Round(ci1_text[1], 2)) * ", " * .(Round(ci1_text[2], 2)) * "]      "
      )
    )),
    col = "black",
    cex = text.cex,
    adj = c(0, 1)
  )

  text(
    x = usr[1] + x.adjust,
    y = y_positions[3],
    labels = as.expression(bquote(
      atop(
        "Between arms using " * .(labels[[2]]) * ":",
        "[" * .(Round(ci2_text[1], 2)) * ", " * .(Round(ci2_text[2], 2)) * "]"
      )
    )),
    col = "black",
    cex = text.cex,
    adj = c(0, 1)
  )

  text(
    x = usr[1] + x.adjust,
    y = y_positions[4],
    labels = as.expression(bquote(
      atop(
        "Between " * .(labels[[1]]) * " and " * .(labels[[2]]) * ":",
        "[" * .(Round(ci_between_text[1], 2)) * ", " * .(Round(ci_between_text[2], 2)) * "]"
      )
    )),
    col = "black",
    cex = text.cex,
    adj = c(0, 1)
  )

  mtext(main, side = 3, line = 0.51, adj = 0, cex = 1.8)
}



#' Plot a single normal distribution
#'
#' @param mean Mean of the distribution.
#' @param sd Standard deviation of the distribution.
#' @param ylab Y-axis label.
#' @param main Plot title.
plot_one_normal_distribution <- function(mean, sd, ylab = "", main){

  # Compute the 95% credible interval for the normal distribution
  lower <- mean - 1.96 * sd
  upper <- mean + 1.96 * sd

  # Create a sequence of x values around mu (covering ±4 sigma)
  x <- seq(mean - 4 * sd, mean + 4 * sd, length.out = 200)
  y <- dnorm(x, mean = mean, sd = sd)
    plot(x, y, type = "l",
         xlim = range(x),    # from min(x) to max(x)
         xaxs = "i",         # remove default padding
         xaxt = "n",         # no automatic x-axis
         xlab = "Difference in ATEs",
         ylab = ylab)


  tickVals <- pretty(range(x))
  tickVals <- sort(unique(c(tickVals, range(x))))

  axis(1, at = tickVals, labels = sprintf("%.1f", tickVals))

  # --- Use segments() to draw dashed lines only up to the curve ---
  y_lower <- dnorm(lower, mean = mean, sd = sd)
  y_upper <- dnorm(upper, mean = mean, sd = sd)

  segments(x0 = lower, y0 = 0, x1 = lower, y1 = y_lower,
           col = scales::alpha("red", alpha = 0.5), lty = 2, lwd = 3)
  segments(x0 = upper, y0 = 0, x1 = upper, y1 = y_upper,
           col = scales::alpha("red", alpha = 0.4), lty = 2, lwd = 3)

  # Shade the area within the credible interval
  idx <- x >= lower & x <= upper
  polygon(c(lower, x[idx], upper), c(0, y[idx], 0),
          col = scales::alpha("red", alpha = 0.1), border = NA)

  # Annotate the credible interval endpoints
  text(lower, 1e-1, labels = sprintf("%.2f", lower),
       pos = 4, col = "red")
  text(upper, 1e-1, labels = sprintf("%.2f", upper),
       pos = 2, col = "red")
  mtext(main, side = 3, line = 0.51, adj = 0, cex = 1.8)
}


# solve quantile
#' Compute selected quantiles
#'
#' @param x Numeric vector.
quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}


# log likelihood of the JSTVC in the first-level
#' Compute the log-likelihood
#'
#' @param data Input data object.
#' @param Para.List Model parameter list.
#' @param yts.minus.xts Response residual matrix after fixed effects.
#' @param slope.Hv.Zg Spatiotemporal fitted component.
#' @param slope.Au.Rg Additional fitted component.
loglik <- function(data, Para.List, yts.minus.xts, slope.Hv.Zg, slope.Au.Rg = 0)
{

  Nt <- data[[1]]$Nt
  n <- 0
  Py <- length(names(data))
  for(py in 1:Py){
    n <- n + data[[py]]$n
  }

  mu.sigma.sq <- vector()
  for (py in 1:Py){
    mu.sigma.sq[data[[py]]$index.y] <- Para.List[[py]]$obs.sigma.sq$mu.sigma.sq
  }

  inv.mu.sigma.sq <- spam::diag.spam(1/mu.sigma.sq)
  mu.sigma.sq <-  spam::diag.spam(mu.sigma.sq)

  Y_ts <- yts.minus.xts - slope.Hv.Zg - slope.Au.Rg

  # ylikelihood
  Loglik.y.1 <- spam::determinant.spam(mu.sigma.sq)$modulus*n*Nt


  Loglik.y.2 <- sapply(seq_len(Nt), function(t){
    t(Y_ts[t, ]) %*% inv.mu.sigma.sq %*% Y_ts[t, ]
  }, simplify = "matrix") %>% sum()

  Loglik <- (- Loglik.y.1 - Loglik.y.2)/2
  return(as.vector(Loglik))
}


#' Evaluate the fixed-effect contribution
#'
#' @param tsv.x Covariate array or list.
#' @param alpha Coefficient vector.
#' @param self Logical; whether to use the self-covariate layout.
Xts.beta.fun <- function(tsv.x, alpha, self = FALSE){
  if(!self){
    Nt <- length(tsv.x)
    if(length(as.vector(alpha)) > 1){
      temp <- sapply(seq_len(Nt), function(t)
        matrix(as.vector(alpha),
               nrow = 1) %*% tsv.x[[t]]
        , simplify = "matrix") #%>% t()
    }else{
      temp <- sapply(seq_len(Nt), function(t)
        tsv.x[[t]]*as.vector(alpha)
        , simplify = "matrix")

    }
    x_ts <- matrix(NA, nrow = Nt, ncol = ncol(tsv.x[[1]]))
    for(t in 1:Nt){
      x_ts[t, ] <- as.vector(temp[[t]]);
    }
  }else{
    Nt <- dim(tsv.x)[[2]]
    if(length(as.vector(alpha)) > 1){
      x_ts <- sapply(seq_len(Nt), function(t)
        matrix(as.vector(alpha),
               nrow = 1) %*% tsv.x[, t, ]
        , simplify = "matrix") %>% t()
    }else{
      x_ts <- sapply(seq_len(Nt), function(t)
        tsv.x[, t, ]*as.vector(alpha)
        , simplify = "matrix") %>% t()

    }
  }
  return(x_ts)
}



#' Transform data lists into stacked matrices
#'
#' @param data Input data object.
#' @param E_inv.sigma.sq Inverse observation variance weights.
#' @param E_inv.obs.delta Additional observation weights.
tranform.list.to.matrix <- function(data, E_inv.sigma.sq = 1, E_inv.obs.delta = 1){


  name.y <- names(data)
  Py <- length(name.y)
  Nt <- data[[1]]$Nt
  n <- 0
  for(py in 1:Py){
    n <- n + data[[py]]$n
  }

  tsv.y <- matrix(NA, nrow = Nt, ncol = n)


  tsv.x <- tsv.sx <- list()
  if(length(E_inv.obs.delta) < (Nt*Py)){
    E_inv.obs.delta <- matrix(E_inv.obs.delta, nrow = Nt, ncol = Py)
  }
  E_inv.sigma.sq <- E_inv.sigma.sq*E_inv.obs.delta

  for(t in 1:Nt) {
    y <- NULL
    for(py in 1:Py){
      y <- c(y, (data[[py]]$Y_ts[t, ]))
    }
    colnames(tsv.y) <- names(y)
    tsv.y[t, ] <- as.vector(y)


    temp.x <- tem.X_ts <- NULL
    for(py in 1:Py){
      if(!is.null(data[[py]]$X_ts)){
        if(data[[py]]$Px == 1){
          tem.X_ts <- matrix(data[[py]]$X_ts[, t, ], nrow = 1, ncol = data[[py]]$n)
        }else{
          tem.X_ts <- data[[py]]$X_ts[, t, ]
        }

      }
      temp.x <- cbind(temp.x, E_inv.sigma.sq[t, py] * tem.X_ts)
    }


    if(!is.null(temp.x)){
      tsv.x[[t]] <- spam::as.spam(as.matrix(temp.x))
    }


  }
  if(length(tsv.x) == 0){
    tsv.x <- NULL
  }
  return(list(tsv.y = tsv.y, tsv.x = tsv.x))
}


#' Initialize EnKF inputs
#'
#' @param data Input training data object.
#' @param Para.List Model parameter list.
#' @param test Optional test data object.
#' @param IDE Logical; whether to use IDE transitions.
#' @param Hv.Zg Current fitted spatiotemporal effect.
#' @param Au.Rg Current fitted additional effect.
IniEnKF <- function (data, Para.List, test = NULL, IDE = TRUE, Hv.Zg = 0, Au.Rg = 0){
  #----------------------------------------------------------------------------
  #-- Stack Zts and H
  train.data.trans.tsv <- tranform.list.to.matrix(data)
  train.Z.tsv          <- test.Z.tsv <- train.sZ.tsv <- test.sZ.tsv <- NULL
  if(!is.null(test)){
    test.data.trans.tsv <- tranform.list.to.matrix(test)
  }
  if(!is.null(train.data.trans.tsv$tsv.z)){
    train.Z.tsv <- train.data.trans.tsv$tsv.z
    if(!is.null(test)){
      test.Z.tsv  <- test.data.trans.tsv$tsv.z
    }
  }
  if(!is.null(train.data.trans.tsv$tsv.sz)){
    train.sZ.tsv <- train.data.trans.tsv$tsv.sz
    if(!is.null(test)){
      test.sZ.tsv  <- test.data.trans.tsv$tsv.sz
    }
  }

  train.H.basis <- list()
  train.H.basis[[1]] <- spam::as.spam(as.matrix(data[[1]]$Hs))
  if(length(data) > 1){
    for(py in 2:length(data)){
      train.H.basis[[py]] <- spam::as.spam(as.matrix(data[[py]]$Hs))
    }
  }

  test.H.basis <- list()
  if(!is.null(test)){
    test.H.basis[[1]] <- spam::as.spam(as.matrix(test[[1]]$Hs))
    if(length(data) > 1){
      for(py in 2:length(data)){
        test.H.basis[[py]] <- spam::as.spam(as.matrix(test[[py]]$Hs))
      }
    }
  }

  # -- transform process
  Mv <- list()
  # term <- c("Intercept", names(Para.List[[1]]$st.RF))
  if(!is.null(data[[1]]$G_ts)){
    term <- c(dimnames(data[[1]]$G_ts)[[1]])
    for(Pg in 1:(length(term))){
      Mv[[term[Pg]]] <- Para.List[[1]]$st.RF[[term[Pg]]]$rho.v$mu.rho.v
      if(IDE){
        Mv[[term[Pg]]] <- list()
        for(g in 1:data$Grid.infor$summary$res){
          Mv[[term[Pg]]][[g]] <-  as(Para.List$st.RF[[term[Pg]]]$rho.v$mu.rho.v[g] *
                                       fields::Wendland(data$Grid.infor$level[[g]]$BAUs.Dist
                                                        , theta = Para.List$st.RF[[term[Pg]]]$Phi.v$mu.Phi.v[g]
                                                        , dimension = 1, k = 1), "dgCMatrix")


          lamda <- rARPACK::eigs_sym(Mv[[term[Pg]]][[g]], 1)$values
          if(lamda >= 1){
            if(is.complex(lamda[1]))
            {
              Lamda <- max(Mod(as.vector(lamda)))
              cat("\n .............................. \n")
              cat("Max-Eigen of the Mv:", Lamda)
              # print(t2 - t1)
              cat("\n .............................. \n")
              if(Lamda > 1) Mv[[term[Pg]]][[g]] <- Mv[[term[Pg]]][[g]]/(abs(Lamda) + 1E-3)
            }else{
              Lamda <- max(lamda)
              cat("\n .............................. \n")
              cat("Max-Eigen of the Mv:", Lamda)
              # print(t2 - t1)
              cat("\n .............................. \n")
              if(Lamda > 1) Mv[[term[Pg]]][[g]] <- Mv[[term[Pg]]][[g]]/(abs(Lamda) + 1E-3)
            }
          }
        }
      }
    }
  }

  for(py in 1:length(names(data))){
    if(!is.null(data[[py]]$sG_ts)){
      term <- dimnames(data[[py]]$sG_ts)[[1]]
      for(sPg in 1:(length(term))){
        Mv[[paste0(names(data)[py], ".", term[sPg])]] <- Para.List[[py]]$st.sRF[[paste0(names(data)[py], ".", term[sPg])]]$rho.v$mu.rho.v
      }
    }
  }


  yts.minus.xts <- train.data.trans.tsv$tsv.y - Hv.Zg - Au.Rg
  if(!is.null(data[[1]]$X_ts)){
    yts.minus.xts <- train.data.trans.tsv$tsv.y - Xts.beta.fun(
      tsv.x = train.data.trans.tsv$tsv.x,
      alpha = Para.List[[1]]$beta$mu.beta,
      self  = FALSE)


  }

  return(list(yts.minus.xts = yts.minus.xts,
              train.H.basis = train.H.basis,
              test.H.basis = test.H.basis,
              Mv = Mv
  ))
}

#' Construct the observation mapping matrix
#'
#' @param data Input data object.
#' @param test Optional test data object.
make.G.Mat <- function (data = NULL, test = NULL){
  data.G.Mat <- test.G.Mat <- Y.index <- Res.indic <- list()
  if(!is.null(data)){
    Nt     <- data[[1]]$Nt
    n      <- 0
    Py     <- length(names(data))
    nKnots <- vector()
    for(py in 1:Py){
      nKnots[py] <- data[[py]]$nKnots
      n <- n + data[[py]]$n
    }
    sPG  <- 0
    ID   <- NULL
    time <- dimnames(data[[py]]$Y_ts)[[1]]
    for(py in 1:Py){
      sPG <- sPG + data[[py]]$sPg
      ID  <- c(ID, paste0("Y", py, ".", dimnames(data[[py]]$Y_ts)[[2]]))
    }

    if(data[[1]]$Pg > 0){
      pb <- progress::progress_bar$new(format = "pub.G.Mat (train): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = data[[1]]$Pg
                                       , clear = FALSE)
      for(Pg in 1:data[[1]]$Pg){
        pb$tick()
        Sys.sleep(1 / (2*data[[1]]$Pg))
        Y.index[[dimnames(data[[1]]$G_ts)[[1]][Pg]]] <- 1:n
        Res.indic[[dimnames(data[[1]]$G_ts)[[1]][Pg]]] <- 0
        data.G.Mat[[dimnames(data[[1]]$G_ts)[[1]][Pg]]] <- array(0,
                                                                 dim =  c(Nt, n, nKnots),
                                                                 dimnames = list(time, ID,
                                                                                 paste0("Knots.", 1:nKnots)))
        for(py in 1:Py){
          for(t in 1:Nt){
            # G_ts <- NULL
            # for(py in 1:Py){
            # G_ts <- c(G_ts, as.vector(data[[py]]$G_ts[Pg, t, ]))
            G_ts <- c(as.vector(data[[py]]$G_ts[Pg, t, ]))
            data.G.Mat[[dimnames(data[[py]]$G_ts)[[1]][Pg]]][t, data[[py]]$index.y,] <- (G_ts*as.matrix(data[[names(data)[py]]]$Hs))
          }

        }
      }
    }
    if(sPG > 0){
      pb <- progress::progress_bar$new(format  = "self.G.Mat (train): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = sPG
                                       , clear = FALSE)
      for(py in 1:Py){
        if(data[[py]]$sPg > 0){
          for(sPg in 1:data[[py]]$sPg){
            Res.indic[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]] <- py
            Y.index[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]] <- data[[py]]$index.y
            pb$tick()
            Sys.sleep(1 / (2*sPG))
            data.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]] <-
              array(0, dim   =  c(Nt, n, nKnots[py]), dimnames = list(time, ID, paste0("Knots.", 1:nKnots[py])))
            for(t in 1:Nt){
              # G_ts                     <- rep(0, n)
              # G_ts[data[[py]]$index.y] <- as.vector(data[[py]]$sG_ts[sPg, t, ])
              G_ts                     <- as.vector(data[[py]]$sG_ts[sPg, t, ])
              data.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]][t, data[[py]]$index.y,] <-
                (G_ts*as.matrix(data[[names(data)[py]]]$Hs))
            }
          }
        }
      }
    }
  }else{
    data.G.Mat <- Y.index <- Res.indic <- NULL
  }
  if(!is.null(test)){
    Nt     <- test[[1]]$Nt
    # nKnots <- test[[1]]$nKnots
    n.test <- 0
    Py     <- length(names(test))
    for(py in 1:Py){
      n.test <- n.test + test[[py]]$n
    }

    time.test <- ID.test <- NULL
    time.test <- c(dimnames(test[[1]]$Y_ts)[[1]])
    sPG       <- 0
    for(py in 1:Py){
      sPG <- sPG + test[[py]]$sPg
      ID.test <- c(ID.test, dimnames(test[[py]]$Y_ts)[[2]])
    }
    if(test[[1]]$Pg > 0){
      pb <- progress::progress_bar$new(format = "pub.G.Mat (test): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = test[[1]]$Pg
                                       , clear = FALSE)

      for(Pg in 1:test[[1]]$Pg){
        pb$tick()
        Sys.sleep(1 / (2*test[[1]]$Pg))
        test.G.Mat[[dimnames(test[[1]]$G_ts)[[1]][Pg]]] <- array(0,
                                                                 dim =  c(Nt, n.test, nKnots),
                                                                 dimnames = list(time.test, ID.test,
                                                                                 paste0("Knots.", 1:nKnots)))
        for(py in 1:Py){
          for(t in 1:Nt){
            #  G_ts <- NULL
            # G_ts <- c(G_ts, as.vector(test[[py]]$G_ts[Pg, t, ]))
            # for(py in 1:Py){
            G_ts <- c(as.vector(test[[names(test)[py]]]$G_ts[Pg, t, ]))
            test.G.Mat[[dimnames(test[[1]]$G_ts)[[1]][Pg]]][t, test[[names(test)[py]]]$index.y, ] <- G_ts*
              as.matrix(test[[names(test)[py]]]$Hs)
          }

        }
      }
    }
    if(sPG > 0){
      pb <- progress::progress_bar$new(format = "self.G.Mat (test): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = sPG
                                       , clear = FALSE)
      for(py in 1:Py){
        if(test[[py]]$sPg > 0){
          for(Pg in 1:test[[py]]$sPg){
            pb$tick()
            Sys.sleep(1 / (2*sPG))
            test.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][Pg])]] <- array(0,
                                                                                                     dim =  c(Nt, n.test, nKnots[py]),
                                                                                                     dimnames = list(time.test, ID.test,
                                                                                                                     paste0("Knots.", 1:nKnots[py])))
            for(t in 1:Nt){
              # G_ts                     <- rep(0, n.test)
              # G_ts[test[[py]]$index.y] <- as.vector(test[[py]]$sG_ts[Pg, t, ])
              G_ts <- as.vector(test[[py]]$sG_ts[Pg, t, ])
              test.G.Mat[[paste0(names(test)[py], ".", dimnames(test[[names(test)[py]]]$sG_ts)[[1]][Pg])]][t,test[[names(test)[py]]]$index.y,] <- G_ts*as.matrix(test[[names(test)[py]]]$Hs) #ini.IEnKF.Input$test.H.basis
            }
          }
        }

      }
    }

  }else{
    test.G.Mat <- NULL
  }
  return(list(data.G.Mat = data.G.Mat, test.G.Mat = test.G.Mat, Y.index = Y.index, Res.indic = Res.indic))
}


#' Compute an empirically tapered covariance matrix
#'
#' @param x Matrix of ensemble draws for the first field.
#' @param y Matrix of ensemble draws for the second field.
#' @param taper Tapering matrix used both as a sparsity mask and weights.
#' @param sym Logical; whether to return a symmetric matrix.
#' @param nThreads Number of threads passed to the C++ helper.
Empical_Cov_taper <- function(x, y, taper, sym = TRUE, nThreads = 1) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  taper <- as.matrix(taper)

  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of ensemble draws.")
  }
  if (nrow(x) != nrow(taper) || nrow(y) != ncol(taper)) {
    stop("Dimensions of `taper` must match `nrow(x)` by `nrow(y)`.")
  }

  idx <- which(taper != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    out <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    if (sym && nrow(out) == ncol(out)) {
      out <- (out + t(out)) / 2
    }
    return(out)
  }

  cov_vec <- Empical_Cov_taper_c(
    x_r = x,
    y_r = y,
    mu_x_r = rowMeans(x),
    mu_y_r = rowMeans(y),
    Index_r = matrix(as.integer(t(idx)), ncol = 2),
    n_r = as.integer(nrow(x)),
    m_r = as.integer(nrow(y)),
    E_r = as.integer(ncol(x)),
    l_r = as.integer(nrow(idx)),
    nThreads_r = as.integer(nThreads)
  )$R

  out <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  out[idx] <- as.vector(cov_vec) * taper[idx]

  if (sym && nrow(out) == ncol(out)) {
    out <- (out + t(out)) / 2
  }

  out
}


#' Build state covariance summaries
#'
#' @param data Input data object.
#' @param Ks State ensemble summaries.
#' @param spTaper Spatial tapering object.
#' @param cores Number of CPU cores to use.
#' @param py Response index.
sfunc <- function(data, Ks, spTaper = NULL, cores = 1, py = 1){
  Pt_1_1 =  parallel::mcmapply(seq_len(data[[py]]$Nt + 1),
                               FUN = function(t){

                                s11 <- spam::as.spam(spam::tcrossprod.spam(Ks$Vt.mean[t,
                                                                                      data[[py]]$Grid.infor$summary$g.index[[1]]])) + spam::as.spam(
                                                                                        cov(t(Ks$Vt.ens[t, data[[py]]$Grid.infor$summary$g.index[[1]], ]),
                                                                                            t(Ks$Vt.ens[t, data[[py]]$Grid.infor$summary$g.index[[1]], ]))*
                                                                                          Ks$rho.time[paste(0)] *spTaper$taper[[1]]);

                                 if(data[[py]]$Grid.infor$summary$res > 1){
                                   for(g in 2:data[[py]]$Grid.infor$summary$res){
                                     s11 <- spam::bdiag.spam(s11,
                                                             spam::tcrossprod.spam(Ks$Vt.mean[t,
                                                                                              data[[py]]$Grid.infor$summary$g.index[[g]]]) + spam::as.spam(
                                                                                                Empical_Cov_taper(Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                                                        Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                                                        Ks$rho.time[paste(0)] *spTaper$taper[[g]],
                                                                                                                        sym = T, nThreads = cores))
                                     );
                                   }
                                 };
                                 return(s11)},
                               SIMPLIFY = F, mc.cores = 1)


  St_1_0 =  parallel::mcmapply(seq_len(data[[py]]$Nt),
                               FUN = function(t){
                                 s10 <- spam::as.spam(spam::tcrossprod.spam(Ks$Vt.mean[t + 1,data[[py]]$Grid.infor$summary$g.index[[1]]],
                                                                            Ks$Vt.mean[t,data[[py]]$Grid.infor$summary$g.index[[1]]])) +
                                   spam::as.spam(
                                     as.matrix(Empical_Cov_taper(Ks$Vt.ens[t + 1,data[[py]]$Grid.infor$summary$g.index[[1]], ],
                                                                      Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[1]], ],
                                                                      Ks$rho.time[paste(1)] *as.matrix(spTaper$taper[[1]]),
                                                                      sym = F,
                                                                      nThreads = cores)));
                                 if(data[[py]]$Grid.infor$summary$res > 1){
                                   for(g in 2:data[[py]]$Grid.infor$summary$res){
                                     s10 <- spam::bdiag.spam(s10,
                                                             spam::tcrossprod.spam(Ks$Vt.mean[t + 1,data[[py]]$Grid.infor$summary$g.index[[g]]],
                                                                                   Ks$Vt.mean[t,data[[py]]$Grid.infor$summary$g.index[[g]]]) +
                                                               spam::as.spam(
                                                                 as.matrix(Empical_Cov_taper(Ks$Vt.ens[t + 1,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                         Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                         Ks$rho.time[paste(1)] *as.matrix(spTaper$taper[[g]]),
                                                                                         sym = F,
                                                                                         nThreads = cores)))
                                     );
                                   }
                                 };
                                 return(s10)},
                               SIMPLIFY = F, mc.cores = 1)
  temp <- 0
  for(t in (2:data[[py]]$Nt)){
    temp <- temp + Pt_1_1[[t]]
  }

  S00 <- temp + Pt_1_1[[1]]
  S11 <- temp + Pt_1_1[[data[[py]]$Nt + 1]]
  S10 <- 0
  for(t in seq_len(data[[py]]$Nt)){
    S10 <- S10 + St_1_0[[t]]#
  }
  rm(St_1_0, temp);gc()




  return(S = list(S00 = S00, S11 = S11, S10 = S10,
                  V.t00 = Pt_1_1[[1]]
  ))
}



#' Objective for variance-parameter optimization
#'
#' @param X Optimization parameter vector.
non_f_var_alpha <- function(X){
  Adj.Mat <- exp(-(G/mu.Phi.v[1])) #

  mat.x <- as.matrix(var.covariable)

  mu.var <- exp(mat.x %*% c(X))
  Adj.Mat <- mu.tau.sq*solve(diag(as.vector(mu.var^(0.5))) %*% Adj.Mat %*% diag(as.vector(mu.var^(0.5))))

  log.det.D <- Matrix::determinant(Adj.Mat)$modulus/2

  f     <-
    (sub.Nt + 1) *log.det.D - (sum(spam::diag.spam(Adj.Mat %*% S11))  +
                                      sum(spam::diag.spam(Adj.Mat %*% V.t00)) -
                                      2*Rho*sum(spam::diag.spam(Adj.Mat %*% S10)) +
                                      Rho.sq*sum(spam::diag.spam(Adj.Mat %*% S00)))/2 -
    0.5*matrix(X - LR.coef.prior$mu.delta, nrow = 1) %*% solve(LR.coef.prior$sigma.sq) %*%
    matrix(X - LR.coef.prior$mu.delta, ncol = 1)
  rm(Adj.Mat);
  return(-f)
}


#' Objective for tau-square KL optimization
#'
#' @param X Optimization parameter vector.
non_f_tau.sq.kl <- function(X){
  a <- exp(X[1]) + 2.1
  b <- exp(X[2])
  E_tau_sq <- b / (a - 1)

  term.1 <- - 0.5*(Upp - Low)* QS.Sum / (1 + E_tau_sq) + 0.5* sub.n*sub.Nt*(log(Low + (Upp - Low)/(1 + E_tau_sq)))
  term.2 <- -(prior.tau.a - 1) * (log(b) - digamma(a)) - prior.tau.b * a / b
  term.3 <- (a - 1) * (log(b) - digamma(a)) - a * log(b) + lgamma(a) + a
  elbo <- term.1 + term.2 + term.3
  return(-elbo)
}


#' Objective for range optimization
#'
#' @param X Optimization parameter vector.
non_f_range <- function(X){
  Hs.fun <- function(X, residuals, Hdist, sigma.sq, Vt.mean, index.y){
    phi <- exp(X[1])
    Hs.temp <- exp(-1*Hdist/phi)
    gg <- Hs.temp/apply(Hs.temp, 1, sum)
    Hs <- as(gg, "sparseMatrix")
    Yt <-  as.vector(residuals[, index.y] - t(Hs %*% t(Vt.mean)))
    return(-Yt%*%Yt/(2*sigma.sq))
  }

  phi <- exp(X[1])
  Adj.Mat <- exp(-G/phi)

  mat.x <- as.matrix(var.covariable)
  mu.var <- exp(mat.x %*% mu.alpha)
  Adj.Mat <- mu.tau.sq*solve(diag(as.vector(mu.var^(0.5))) %*% Adj.Mat %*% diag(as.vector(mu.var^(0.5))))
  log.det.D <- Matrix::determinant(Adj.Mat)$modulus/2

  f     <- Hs.fun(X, residuals, Hdist, sigma.sq, Vt.mean, index.y) -
    (sub.Nt + 1) *log.det.D - (sum(spam::diag.spam(Adj.Mat %*% S11))  +
                                      sum(spam::diag.spam(Adj.Mat %*% V.t00)) -
                                      2*Rho*sum(spam::diag.spam(Adj.Mat %*% S10)) +
                                      Rho.sq*sum(spam::diag.spam(Adj.Mat %*% S00)))/2
  rm(Adj.Mat);
  return(-f)
}

###################################################################
#' Validate spatiotemporal predictions
#'
#' @param z Observed values.
#' @param zhat Predicted mean values.
#' @param sigma Prediction uncertainty values.
#' @param zhat.Ens Predictive ensemble samples.
#' @param names Logical; whether to print metric names.
#' @param CC Logical; whether to compute the full metric set.
spT_validation <- function (z = NULL, zhat = NULL, sigma = NULL,
                            zhat.Ens = NULL, names = FALSE, CC = TRUE){
  VMSE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), 4)
  }
  ## root mean square error
  RMSE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), 4)
  }
  MAE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  MAPE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))/abs(z)
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u) * 100, 4)
  }
  ## normalised mean gross error
  NMGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    NMGE <- abs(c(x[[2]] - x[[1]]))/sum(x[[1]])
    round(mean(NMGE), 4)
  }

  BIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  rBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u) * mean(z, na.rm = TRUE)), 4)
  }
  ## normalised mean bias
  nBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    round(sum(x, na.rm = TRUE) * 100 / sum(z, na.rm = TRUE), 4)
  }

  rMSEP <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat, na.rm = TRUE) - z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
  }
  ## correlation coefficient
  Coef <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    Coef <- suppressWarnings(cor(x[[2]], x[[1]])) ## when SD=0; will return NA
    round(Coef, 4)
  }

  ##  Coefficient of Efficiency
  COE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    COE <- 1 - sum(abs(x[[2]] - x[[1]])) / sum(abs(x[[1]] - mean(x[[1]])))
    round(COE, 4)
  }
  ## fraction within a factor of two
  FAC2 <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    ratio <- x[[2]]/x[[1]]
    len <- length(ratio)
    if (len > 0) {
      FAC2 <- length(which(ratio >= 0.5 & ratio <= 2)) / len
    } else {
      FAC2 <- NA
    }
    round(FAC2, 4)
  }
  ## mean gross error
  MGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    MGE <- round(mean(abs(x[[2]] - x[[1]])), 4)
  }

  ##  Index of Agreement
  IOA <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)

    LHS <- sum(abs(x[[2]] - x[[1]]))
    RHS <- 2 * sum(abs(x[[1]] - mean(x[[1]])))

    if (LHS <= RHS) IOA <- 1 - LHS / RHS else IOA <- RHS / LHS - 1
    round(IOA, 4)
  }


  CRPS <<- function(z, zhat, sigma = NA){
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    if(is.na(sigma[1])){
      da <- data.frame(z = as.vector(z), zhat = as.vector(zhat))
      x <- na.omit(da)
      x$sigma <- sd(x$zhat)
    }else{
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat),
                       sigma = as.vector(sigma))
      x <- na.omit(da)
    }
    CRPS <- round(mean(verification::crps(x$z, cbind(x$zhat, x$sigma))$CRPS), 4)
  }


  Interval.score <<- function(z, zhat, sigma = NA, alpha=0.05) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    if(is.na(sigma[1])){
      da <- data.frame(z = as.vector(z), zhat = as.vector(zhat))
      x <- na.omit(da)
      x$sigma <- sd(x$zhat)
    }else{
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat),
                       sigma = as.vector(sigma))
      x <- na.omit(da)
    }

    ## Compute the Interval score:
    hw <- -qnorm(alpha/2) * x$sigma
    scores <- 2 * hw + (2/alpha) * (((x$z - hw) - x$zhat) * (x$zhat < x$z - hw) +
                                      (x$zhat - (x$z + hw)) * (x$zhat > x$z + hw))
    # attr(scores, "score") <- mean(scores, na.rm=TRUE)
    scores <- mean(scores, na.rm=TRUE)
  }

  CRPS.ES <- function(z, zhat, zhat.Ens) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))

    Nt <- nrow(z)
    n <- ncol(z)
    # cat("Nt = ", Nt)
    if(!is.null(zhat.Ens))
    {
      # z <- as.matrix(z)
      test.index <- which(is.na(z), arr.ind = T)
      if(nrow(test.index) >0){ z[test.index] <- 0 }
      CRPS <- ES <- matrix(NA, nrow = nrow(z), ncol = n)
      for(t in 1:Nt)
      {
        CRPS[t, ] <- crps_sample(z[t, ], zhat.Ens[t,,])
        ES[t, ] <- es_sample(z[t, ], zhat.Ens[t,,])
      }
      CRPS[test.index] <-  ES[test.index] <- NA
      #crps
      CRPS <- round(mean(CRPS, na.rm = TRUE), 4)
      ES <- round(mean(ES, na.rm = TRUE), 4)
    }else{
      CRPS <- ES <- vector()
      for(t in 1:Nt)
      {
        test.index <- which(is.na(z[t, ]))%>% as.vector()
        if(length(test.index) == 0 ) {
          z0 <- z[t, ] %>% as.vector()
          zt <- zhat[t, ] %>% as.matrix()
        }else{
          z0 <- z[t, -test.index] %>% as.vector()
          zt <- zhat[t, -test.index] %>% as.matrix()
        }
        CRPS[t] <- mean(crps_sample(z0, zt))
        ES[t] <-  es_sample(z0, zt)
      }
      CRPS <- round(mean(CRPS, na.rm = TRUE), 4)
      ES <- round(mean(ES, na.rm = TRUE), 4)
    }
    return(list(CRPS = CRPS, ES = ES))
  }

  CRPS.samples <<- function(true, samples, weights=-1) {
    true <- as.vector(true)
    Ne <- dim(samples)[3]
    if(!is.null(Ne)){
      Sample <- samples
      samples <- NULL
      for(e in 1:Ne){
        samples <- cbind(samples, as.vector(Sample[,, e]))
      }
    }


    m <- length(true)

    if((is.vector(samples) & m > 1) || ncol(samples)==1) {

      expect1 = abs(samples-true)
      expect2 = 0

    }else{

      if(is.vector(samples)){ samples=matrix(samples,nrow=1,ncol=length(samples)) }

      n=ncol(samples)

      if(sum(weights==-1)>0) {
        weights=matrix(1/n,nrow=nrow(samples),ncol=n)
      }

      expect1 = rowSums(abs(samples-true)*weights)

      sortind=t(apply(samples,1,order))

      for(i in 1:m) {
        samples[i,]=samples[i,sortind[i,]]
        weights[i,]=weights[i,sortind[i,]]
      }

      cum.weights=t(apply(weights,1,cumsum))[,1:(n-1)]

      if(m==1) {
        expect2 = sum(cum.weights*(1-cum.weights)*(samples[2:n]-samples[1:(n-1)]))
      } else {
        expect2 = rowSums(cum.weights*(1-cum.weights)*(samples[,2:n]-samples[,1:(n-1)]))
      }

    }

    return(round(mean(expect1 - expect2), 4))

  }
  Ignorance_score <<- function(z, zhat, zhat.Ens){
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    if(!is.null(zhat.Ens)){
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat),
                       zhat.Ens = as.vector(apply(zhat.Ens,
                                                  c(1, 2), sd)))
      x <- na.omit(da)
      IGN <- verification::crps(x$z,
                                cbind(x$zhat,
                                      x$zhat.Ens))$IGN
    }else{
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat))
      x <- na.omit(da)
      IGN <- verification::crps(x$z,
                                cbind(x$zhat, sd(x$zhat)))$IGN
    }
    return(IGN)
  }
  if (names == TRUE) {
    cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Normalized Mean Error(NME: %) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Normalized Mean Bias(NMB: %) \n Relative Mean Separation (rMSEP) \n Correlation coefficient (Coef.) \n The fraction of predictions within a factor of two of observations (FAC2) \n Continuous ranked probability score (CRPS) based on sample mean and standard deviation \n CRPS based on empirical distribution function \n Energy score (ES) based on empirical distribution function \n##\n")
  }
  if((!is.null(z))&(!is.null(zhat))){

    out <- NULL
    if(CC){
      out$MSE <- VMSE((z), (zhat))
    }

    out$RMSE <- RMSE((z), (zhat))
    if(CC){
      out$MAPE <- MAPE((z), (zhat))
      out$BIAS <- BIAS((z), (zhat))
      out$rBIAS <- rBIAS((z), (zhat))
    }

    out$Coef <- Coef((z), (zhat))
    out$FAC2 <- FAC2(z, zhat)
    out$CRPS <- CRPS(z, zhat, sigma = sigma)
    out$MAE <- MAE((z), (zhat))
    # out$INS <- Interval.score(z, zhat, sigma = sigma)
    if(!is.null(zhat.Ens)){
      out$CRPS.samples <- CRPS.samples(z, zhat.Ens)
    }else{
      out$CRPS.samples <- CRPS.samples(as.vector(z), as.vector(zhat))
    }

    # out$IGN <- Ignorance_score(z, zhat, zhat.Ens)

    if(CC){
      if(!is.null(zhat.Ens))
      {
        R <- CRPS.ES(z, zhat, zhat.Ens)
      }else{
        R <- CRPS.ES(z, zhat, zhat.Ens = NULL)
      }

      out$CRPS.sample <- R$CRPS
      out$ES.sample <- R$ES
      out$NMGE <- NMGE((z), (zhat))
      out$MGE <- MGE(z, zhat)
      out$IOA <- IOA(z, zhat)
      out$NMB <- nBIAS((z), (zhat))
      out$rMSEP <- rMSEP((z), (zhat))
      out$COE <- COE(z, zhat)
    }
    unlist(out)
  }else{
    return(0)
  }
}


######################################################################
#                      plot  scatter
######################################################################
#' Format rounded numeric values
#'
#' @param x Numeric values.
#' @param n Number of decimal places.
Round <- function(x, n = 2)
{
  return(format(round(x, n), nsmall = n))
}

#Summary for competing methods
#' Summarize model performance
#'
#' @param beta_mean Matrix of posterior means.
#' @param beta_var Matrix of posterior variances.
#' @param beta_true Vector of true coefficient values.
summarize_model <- function(beta_mean, beta_var, beta_true) {
  performance <- calc_performance(beta_mean, beta_true)
  coverage_prob <- round(
    colMeans(compute_coverage(beta_mean, beta_var, beta_true)) * 100,
    5
  )

  list(
    coefficient = performance$Details$Coefficient,
    bias = performance$Details$Bias,
    sd = performance$Details$SD,
    mse = performance$Details$MSE,
    ec = coverage_prob
  )
}

#' Collect competing-model results
#'
#' @param save.Tab.from Directory containing saved result files.
collect_competing_results <- function(save.Tab.from) {
  model_results <- list(
    nonEffect = list(beta_mean = NULL, beta_var = NULL),
    subRegion = list(beta_mean = NULL, beta_var = NULL),
    subArm = list(beta_mean = NULL, beta_var = NULL)
  )

  files <- list.files(save.Tab.from)
  files <- files[order(as.numeric(sub("sim_(\\d+)\\.RData", "\\1", files)))]

  for (f in files) {
    load(file.path(save.Tab.from, f))

    model_results$nonEffect$beta_mean <- rbind(
      model_results$nonEffect$beta_mean,
      t(JSTVC.nonEffect.beta)
    )
    model_results$subRegion$beta_mean <- rbind(
      model_results$subRegion$beta_mean,
      t(JSTVC.subRegion.beta)
    )
    model_results$subArm$beta_mean <- rbind(
      model_results$subArm$beta_mean,
      t(JSTVC.subArm.beta)
    )

    model_results$nonEffect$beta_var <- rbind(
      model_results$nonEffect$beta_var,
      t(JSTVC.nonEffect.sd^2)
    )
    model_results$subRegion$beta_var <- rbind(
      model_results$subRegion$beta_var,
      t(JSTVC.subRegion.sd^2)
    )
    model_results$subArm$beta_var <- rbind(
      model_results$subArm$beta_var,
      t(JSTVC.subArm.sd^2)
    )
  }

  model_results
}

#' Create a summary table for competing models
#'
#' @param save.Tab.from Directory containing saved result files.
#' @param beta_true Vector of true coefficient values.
create_competing_summary_table <- function(save.Tab.from,
                                           beta_true = c(5, rep(-1, 10))) {
  model_results <- collect_competing_results(save.Tab.from)

  nonEffect <- summarize_model(
    model_results$nonEffect$beta_mean,
    model_results$nonEffect$beta_var,
    beta_true
  )
  subRegion <- summarize_model(
    model_results$subRegion$beta_mean,
    model_results$subRegion$beta_var,
    beta_true
  )
  subArm <- summarize_model(
    model_results$subArm$beta_mean,
    model_results$subArm$beta_var,
    beta_true
  )

  data.frame(
    beta = nonEffect$coefficient,
    nonEffect.bias = nonEffect$bias,
    subRegion.bias = subRegion$bias,
    subArm.bias = subArm$bias,
    nonEffect.sd = nonEffect$sd,
    subRegion.sd = subRegion$sd,
    subArm.sd = subArm$sd,
    nonEffect.mse = nonEffect$mse,
    subRegion.mse = subRegion$mse,
    subArm.mse = subArm$mse,
    nonEffect.ec = nonEffect$ec,
    subRegion.ec = subRegion$ec,
    subArm.ec = subArm$ec
  )
}

#' Write the competing-model summary table
#'
#' @param save.Tab.from Directory containing saved result files.
#' @param save.Tab.to Output path for the summary workbook.
#' @param beta_true Vector of true coefficient values.
write_competing_summary <- function(save.Tab.from,
                                    save.Tab.to,
                                    beta_true = c(5, rep(-1, 10))) {
  da <- create_competing_summary_table(
    save.Tab.from = save.Tab.from,
    beta_true = beta_true
  )
  writexl::write_xlsx(da, path = save.Tab.to)
  invisible(da)
}

# summary for JSTVC
#' Create a wide comparison table
#'
#' @param Competing.method Data frame of competing-method summaries.
#' @param JSTVC_ie Summary object for the JSTVC without IE component.
#' @param JSTVC Summary object for the JSTVC model.
#' @param metric Metric name to extract.
#' @param models Model labels.
#' @param n_beta Number of coefficient columns.
#' @param digits Number of digits to display.
#' @param Opt Logical; whether to highlight optimal values.
create_wide_table <- function(Competing.method,
                              JSTVC_ie = NULL,
                              JSTVC,
                              metric = "Bias",
                              models = c("JSTVC-xi", "Sub-AR1", "Arm-AR1", "JSTVC_IE", "JSTVC"),
                              n_beta = 11,
                              digits = 3,
                              Opt = TRUE) {
  beta_names <- paste0("Beta", 1:n_beta)

  da_wide <- data.frame(
    Metric = metric,
    Model = models
  )

  cols <- grep(paste0("\\.", tolower(metric), "$"), names(Competing.method), value = TRUE)
  beta_matrix <- NULL

  for (l in seq_along(cols)) {
    x <- Competing.method[[cols[l]]]
    beta_matrix <- rbind(beta_matrix, t(formatC(x, format = "f", digits = digits)))
  }

  if (!is.null(JSTVC_ie)) {
    cols_jstvc_ie <- grep(metric, names(JSTVC_ie$Details), value = TRUE)
    jstvc_ie_matrix <- cbind(
      do.call(rbind, lapply(cols_jstvc_ie, function(col) {
        t(formatC(JSTVC_ie$Details[, col], format = "f", digits = digits))
      })),
      "", ""
    )
    beta_matrix <- rbind(beta_matrix, jstvc_ie_matrix)
  }

  cols_jstvc <- grep(metric, names(JSTVC$Details), value = TRUE)
  jstvc_matrix <- do.call(rbind, lapply(cols_jstvc, function(col) {
    t(formatC(JSTVC$Details[, col], format = "f", digits = digits))
  }))
  beta_matrix <- rbind(beta_matrix, jstvc_matrix)

  colnames(beta_matrix) <- beta_names
  da_wide[, beta_names] <- beta_matrix

  if (Opt) {
    for (col in 3:(n_beta + 2)) {
      numeric_vals <- abs(as.numeric(da_wide[[col]]))
      best_val <- min(numeric_vals, na.rm = TRUE)
      da_wide[[col]] <- ifelse(
        numeric_vals == best_val,
        sprintf("\\textbf{%.4f}", numeric_vals),
        sprintf("%.4f", numeric_vals)
      )
    }
  }

  da_wide
}

#' Summarize one simulation setting
#'
#' @param competing_path Path to the competing-method summary file.
#' @param save.Tab.from1 Directory containing JSTVC result files.
#' @param save.Tab.from2 Directory containing JSTVC IE result files.
#' @param save.Tab.to Output path for the summary workbook.
#' @param beta_true Vector of true coefficient values.
summarize_one_setting <- function(competing_path,
                                  save.Tab.from1,
                                  save.Tab.from2,
                                  save.Tab.to,
                                  beta_true = c(5, rep(-1, 10))) {
  Competing.method <- readxl::read_xlsx(competing_path)

  JSTVC_beta_mean <- JSTVC_beta_var <- NULL
  JSTVC_IE_beta_mean <- JSTVC_IE_beta_var <- NULL

  JSTVC_files <- list.files(save.Tab.from1)
  JSTVC_files <- JSTVC_files[order(as.numeric(sub("sim_(\\d+)\\.RData", "\\1", JSTVC_files)))]
  for (f in JSTVC_files) {
    load(file.path(save.Tab.from1, f))
    JSTVC_beta_mean <- rbind(JSTVC_beta_mean, t(beta))
    JSTVC_beta_var  <- rbind(JSTVC_beta_var, t(sigma.sq))
  }

  JSTVC_IE_files <- list.files(save.Tab.from2)
  JSTVC_IE_files <- JSTVC_IE_files[order(as.numeric(sub("sim_(\\d+)\\.RData", "\\1", JSTVC_IE_files)))]
  for (f in JSTVC_IE_files) {
    load(file.path(save.Tab.from2, f))
    JSTVC_IE_beta_mean <- rbind(JSTVC_IE_beta_mean, t(beta))
    JSTVC_IE_beta_var  <- rbind(JSTVC_IE_beta_var, t(sigma.sq))
  }

  JSTVC    <- calc_performance(JSTVC_beta_mean, beta_true)
  JSTVC_ie <- calc_performance(JSTVC_IE_beta_mean, beta_true[1:9])

  JSTVC$Details$JSTVC.EC <- round(
    colMeans(compute_coverage(JSTVC_beta_mean, JSTVC_beta_var, beta_true)) * 100, 5
  )
  JSTVC_ie$Details$JSTVC.EC <- round(
    colMeans(compute_coverage(JSTVC_IE_beta_mean, JSTVC_IE_beta_var, beta_true[1:9])) * 100, 5
  )

  models <- c(
    "$\\mathrm{JSTVC}_{-\\xi}$",
    "Sub-$\\mathrm{AR}_{1}$",
    "Arm-$\\mathrm{AR}_{1}$",
    "$\\mathrm{JSTVC}_{-\\mathrm{IE}}$",
    "JSTVC"
  )

  Bias <- create_wide_table(Competing.method, JSTVC_ie = JSTVC_ie, JSTVC,
                            metric = "Bias", models = models, digits = 5)
  SD   <- create_wide_table(Competing.method, JSTVC_ie = JSTVC_ie, JSTVC,
                            metric = "SD", models = models, digits = 5)
  MSE  <- create_wide_table(Competing.method, JSTVC_ie = JSTVC_ie, JSTVC,
                            metric = "MSE", models = models, digits = 5)
  EC   <- create_wide_table(Competing.method, JSTVC_ie = JSTVC_ie, JSTVC,
                            metric = "EC", models = models, digits = 2, Opt = FALSE)

  colnames(Bias) <- colnames(SD) <- colnames(MSE) <- colnames(EC) <-
    c("Metric", "Model", "alpha_0",
      paste0("alpha_1", 1:4), paste0("alpha_2", 1:4),
      "gamma_1", "gamma_2")

  EC$alpha_0 <- NULL
  writexl::write_xlsx(
    list(Statistics = rbind(Bias, SD = SD, MSE = MSE), EC = EC),
    path = save.Tab.to
  )

  invisible(list(Bias = Bias, SD = SD, MSE = MSE, EC = EC))
}
