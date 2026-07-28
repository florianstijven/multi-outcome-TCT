#!/usr/bin/env Rscript

simulate_appendix_c <- function(
    j_values = c(5, 10, 20, 50),
    scale_grid = seq(0, 6, length.out = 81),
    alpha = 0.05,
    output_dir = getwd()
) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for the bivariate normal CDF calculation.")
  }

  normal_pdf <- function(x) dnorm(x)

  probit_mean <- function(mu, a, b, c) {
    c + (1 - c) * pnorm(a * (mu - b) / sqrt(1 + a^2))
  }

  bivariate_normal_cdf <- function(x, y, rho) {
    sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
    as.numeric(mvtnorm::pmvnorm(upper = c(x, y), corr = sigma))
  }

  item_covariance_matrix <- function(mu, a, b, c) {
    j <- length(a)
    mean_vec <- vapply(seq_len(j), function(idx) probit_mean(mu, a[idx], b[idx], c[idx]), numeric(1))
    sigma <- matrix(0, nrow = j, ncol = j)

    for (idx in seq_len(j)) {
      sigma[idx, idx] <- mean_vec[idx] * (1 - mean_vec[idx])
    }

    if (j >= 2) {
      for (idx in seq_len(j - 1)) {
        for (idy in (idx + 1):j) {
          rho <- a[idx] * a[idy] / sqrt((1 + a[idx]^2) * (1 + a[idy]^2))
          x <- a[idx] * (mu - b[idx]) / sqrt(1 + a[idx]^2)
          y <- a[idy] * (mu - b[idy]) / sqrt(1 + a[idy]^2)

          e_phi_phi <- bivariate_normal_cdf(x, y, rho)
          e_eta_eta <- c[idx] * c[idy] + c[idx] * (1 - c[idy]) * mean_vec[idy] + c[idy] * (1 - c[idx]) * mean_vec[idx] +
            (1 - c[idx]) * (1 - c[idy]) * e_phi_phi

          sigma[idx, idy] <- e_eta_eta - mean_vec[idx] * mean_vec[idy]
          sigma[idy, idx] <- sigma[idx, idy]
        }
      }
    }

    list(mean = mean_vec, cov = sigma)
  }

  local_effect_vector <- function(a, b, c, mu0 = 0) {
    (1 - c) * normal_pdf(a * (mu0 - b) / sqrt(1 + a^2)) * a / sqrt(1 + a^2)
  }

  local_parameter_vector <- function(scale, direction) {
    scale * direction
  }

  power_from_ncp <- function(ncp, df, alpha) {
    critical_value <- qchisq(1 - alpha, df = df)
    1 - pchisq(critical_value, df = df, ncp = ncp)
  }

  scenario_parameters <- function(j, scenario_name) {
    index <- seq_len(j)
    if (scenario_name == "scenario_1") {
      list(a = rep(1, j), b = rep(0, j), c = rep(0, j))
    } else if (scenario_name == "scenario_2") {
      list(a = 3 * (index / j), b = rep(0, j), c = rep(0, j))
    } else if (scenario_name == "scenario_3") {
      list(a = 3 * (index / j), b =  2 * (index / j) * (-1)^index, c = rep(0.20, j))
    } else if (scenario_name == "scenario_4") {
      list(a = 1 + 2 * (index / j), b = 3 * ((index - 1) / max(1, j - 1) - 1 / 2), c = rep(0, j))
    } else {
      stop("Unknown scenario: ", scenario_name)
    }
  }

  compute_power_curve <- function(j, scenario_name) {
    params <- scenario_parameters(j, scenario_name)
    moments <- item_covariance_matrix(mu = 0, a = params$a, b = params$b, c = params$c)
    sigma0 <- moments$cov
    d_vec <- local_effect_vector(params$a, params$b, params$c)
    ones <- rep(1, j)

    sum_variance <- as.numeric(2 * crossprod(ones, sigma0 %*% ones))
    sigma0_inv <- solve(sigma0)

    out <- data.frame(
      scenario = scenario_name,
      J = j,
      scale = scale_grid,
      power_t_test = NA_real_,
      power_multivariate = NA_real_,
      ncp_t_test = NA_real_,
      ncp_multivariate = NA_real_,
      stringsAsFactors = FALSE
    )

    for (idx in seq_along(scale_grid)) {
      scale <- scale_grid[idx]
      h_vec <- local_parameter_vector(scale, d_vec)

      ncp_t <- as.numeric(crossprod(ones, h_vec))^2 / sum_variance
      ncp_mv <- as.numeric(crossprod(h_vec, sigma0_inv %*% h_vec))

      out$ncp_t_test[idx] <- ncp_t
      out$ncp_multivariate[idx] <- ncp_mv
      out$power_t_test[idx] <- power_from_ncp(ncp_t, df = 1, alpha = alpha)
      out$power_multivariate[idx] <- power_from_ncp(ncp_mv, df = j, alpha = alpha)
    }

    out
  }

  scenarios <- c("scenario_1", "scenario_2", "scenario_3", "scenario_4")
  results <- do.call(
    rbind,
    lapply(j_values, function(j) {
      do.call(rbind, lapply(scenarios, function(scenario_name) compute_power_curve(j, scenario_name)))
    })
  )

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  csv_file <- file.path(output_dir, "appendix_c_power_curves.csv")
  write.csv(results, csv_file, row.names = FALSE)

  png_file <- file.path(output_dir, "appendix_c_power_curves.png")
  grDevices::png(png_file, width = 2200, height = 2000, res = 220)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(scenarios), length(j_values)), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  for (scenario_name in scenarios) {
    for (j in j_values) {
      panel <- results[results$scenario == scenario_name & results$J == j, ]

      plot(
        panel$scale,
        panel$power_multivariate,
        type = "l",
        lwd = 2,
        col = "#1f77b4",
        ylim = c(0, 1),
        xlab = "local scale s",
        ylab = "Power",
        main = paste(scenario_name, "J =", j)
      )
      lines(panel$scale, panel$power_t_test, lwd = 2, col = "#d62728", lty = 2)
      abline(h = alpha, col = "grey70", lty = 3)
      legend(
        "bottomright",
        legend = c("Item-level chi-squared", "Total-score t-test"),
        col = c("#1f77b4", "#d62728"),
        lty = c(1, 2),
        lwd = 2,
        bty = "n",
        cex = 0.75
      )
    }
  }

  mtext("Appendix C power curves", outer = TRUE, cex = 1.2)

  invisible(list(results = results, csv_file = csv_file, png_file = png_file))
}

if (sys.nframe() == 0) {
  simulate_appendix_c()
}

