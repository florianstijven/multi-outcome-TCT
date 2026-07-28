#!/usr/bin/env Rscript

simulate_local_power_working_models <- function(
    j = 5,
    k = 5,
    h_grid = seq(0, 1, length.out = 80),
    alpha = 0.05
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  time_index <- 0:k
  time_corr <- outer(time_index, time_index, function(i, j) 0.5^abs(i - j))
  outcome_corr <- outer(seq_len(j), seq_len(j), function(i, j) 0.5^abs(i - j))

  sigma_effect <- diag(c(1, 1)) %x% outcome_corr %x% time_corr

  base_effect <- rep(c(rep(0, k + 1), time_index), times = j)

  build_omnibus_contrast_single_outcome <- function(k) {
    contrast <- matrix(0, nrow = k, ncol = 2 * (k + 1))
    contrast[, 2:(k + 1)] <- diag(-1, k)
    contrast[, (k + 3):(2 * (k + 1))] <- diag(1, k)
    contrast
  }

  build_omnibus_contrast_multi_outcome <- function(j, k) {
    diag(1, j) %x% build_omnibus_contrast_single_outcome(k)
  }

    build_linear_contrast_multi_outcome <- function(j, k, time_index, Sigma) {
        jacobian <- matrix(0, nrow = j, ncol = 2 * j * (k + 1))
        for (outcome_idx in seq_len(j)) {
            start_col <- (outcome_idx - 1) * 2 * (k + 1) + (k + 2)
            end_col <- start_col + k
            jacobian[outcome_idx, start_col:end_col] <- time_index
        }
        jacobian = t(jacobian)

    A_contrast <- build_omnibus_contrast_multi_outcome(j, k)

    B_contrast <- t(A_contrast) %*% solve(A_contrast %*% Sigma %*% t(A_contrast)) %*% A_contrast %*% jacobian 
    B_contrast <- t(B_contrast)

    return(B_contrast)
    }

    build_linear_contrast_common <- function(j, k, time_index, Sigma) {
        jacobian <- matrix(0, nrow = 1, ncol = 2 * j * (k + 1))
        for (outcome_idx in seq_len(j)) {
            start_col <- (outcome_idx - 1) * 2 * (k + 1) + (k + 2)
            end_col <- start_col + k
            jacobian[1, start_col:end_col] <- time_index
        }
        jacobian = t(jacobian)

    A_contrast <- build_omnibus_contrast_multi_outcome(j, k)

    B_contrast <- t(A_contrast) %*% solve(A_contrast %*% Sigma %*% t(A_contrast)) %*% A_contrast %*% jacobian
    B_contrast <- t(B_contrast)

    return(B_contrast)
    }

  power_from_ncp <- function(ncp, df, alpha) {
    critical_value <- qchisq(1 - alpha, df = df)
    1 - pchisq(critical_value, df = df, ncp = ncp)
  }

  compute_power_curve <- function(B, h_grid, alpha, Sigma, base_effect) {
     df = nrow(B)

    out <- data.frame(
      h = h_grid,
      power = NA_real_,
      ncp = NA_real_,
      df = df,
      stringsAsFactors = FALSE
    )

    for (idx in seq_along(h_grid)) {
      h <- h_grid[idx]
      delta <- h * base_effect
      ncp_value <- as.numeric(
        t(delta) %*% t(B) %*% solve(B %*% Sigma %*% t(B)) %*% B %*% delta
      )
      out$power[idx] <- power_from_ncp(ncp_value, df = df, alpha = alpha)
      out$ncp[idx] <- ncp_value
    }

    out
  }

  B_omnibus <- build_omnibus_contrast_multi_outcome(j, k)
  B_linear <- build_linear_contrast_multi_outcome(j, k, time_index, sigma_effect)
  B_linear_common <- build_linear_contrast_common(j, k, time_index, sigma_effect)

  results <- list(
    omnibus = compute_power_curve(B_omnibus, h_grid = h_grid, alpha = alpha, Sigma = sigma_effect, base_effect = base_effect),
    linear = compute_power_curve(B_linear, h_grid = h_grid, alpha = alpha, Sigma = sigma_effect, base_effect = base_effect),
    true = compute_power_curve(B_linear_common, h_grid = h_grid, alpha = alpha, Sigma = sigma_effect, base_effect = base_effect)
  )

  return(results)
}

plot_local_power_working_models <- function(    
    j_values = c(3, 6, 12, 24),
    k_values = 5,
    h_grid = seq(0, 0.5, length.out = 80),
    alpha = 0.05,
    output_dir = getwd()) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

png_file <- file.path(output_dir, "local_power_working_models.png")
  grDevices::png(png_file, width = 2200, height = 800, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)

    par(mfrow = c(length(k_values), length(j_values)), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  for (k in k_values) {
    for (j in j_values) {
    results <- simulate_local_power_working_models(j = j, k = k, h_grid = h_grid, alpha = alpha)
      plot(
    h_grid,
    results$omnibus$power,
    type = "l",
    lwd = 2,
    col = "#1f77b4",
    ylim = c(0, 1),
    xlab = "Local effect size h",
    ylab = "Power",
    main = paste("J =", j, " K =", k)
  )
  lines(h_grid, results$linear$power, lwd = 2, col = "#d62728", lty = 2)
  lines(h_grid, results$true$power, lwd = 2, col = "#2ca02c", lty = 3)
  abline(h = alpha, col = "grey70", lty = 3)
  legend(
    "bottomright",
    legend = c("Omnibus test", "Linear working-model test", "Oracle test"),
    col = c("#1f77b4", "#d62728", "#2ca02c"),
    lty = c(1, 2, 3),
    lwd = 2,
    bty = "n"
  )
    }
  }

invisible(list(results = results, png_file = png_file))

# Add same plots, but with relative power, using the omnibus test as the reference.
png_file <- file.path(output_dir, "local_power_working_models_relative.png")
  grDevices::png(png_file, width = 2200, height = 800, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)

    par(mfrow = c(length(k_values), length(j_values)), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  for (k in k_values) {
    for (j in j_values) {

    results <- simulate_local_power_working_models(j = j, k = k, h_grid = h_grid, alpha = alpha)
    results$linear$power_relative <- results$linear$power / results$omnibus$power
    results$true$power_relative <- results$true$power / results$omnibus$power
    results$omnibus$power_relative <- results$omnibus$power / results$omnibus$power
      plot(
    h_grid,
    results$omnibus$power_relative,
    type = "l",
    lwd = 2,
    col = "#1f77b4",
    ylim = c(1, 8),
    xlab = "Local effect size h",
    ylab = "Relative Power",
    main = paste("J =", j, " K =", k)
  )
  lines(h_grid, results$linear$power_relative, lwd = 2, col = "#d62728", lty = 2)
  lines(h_grid, results$true$power_relative, lwd = 2, col = "#2ca02c", lty = 3)
  abline(h = alpha, col = "grey70", lty = 3)
  legend(
    "bottomright",
    legend = c("Omnibus test", "Linear working-model test", "Oracle test"),
    col = c("#1f77b4", "#d62728", "#2ca02c"),
    lty = c(1, 2, 3),
    lwd = 2,
    bty = "n"
  )
    }
  }
    }

if (sys.nframe() == 0) {
  plot_local_power_working_models()
}
