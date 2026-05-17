## Public API for the variance of the empirical mean under a birth-death
## process with Brownian traits.

empirical_mean_variance <- function(birth,
                                    death,
                                    sigma2 = 1,
                                    time_start = 0,
                                    time_end,
                                    time_step,
                                    conditioning = "survival",
                                    n_steps = 1200,
                                    tol_sing = 1e-7,
                                    li2_tol = 1e-12) {
  if (!identical(conditioning, "survival")) {
    stop("Only conditioning = 'survival' is implemented by the original script.", call. = FALSE)
  }

  grid <- .make_public_grid(time_start, time_end, time_step)
  birth_fun <- .as_rate_function(birth, "birth")
  death_fun <- .as_rate_function(death, "death")
  env <- .bd_env_meanvar()

  res <- env$compute_theory_BD_BM_meanvar(
    lambda_fun = birth_fun,
    mu_fun = death_fun,
    sigma2 = sigma2,
    grid = grid,
    n_steps = n_steps,
    tol_sing = tol_sing,
    li2_tol = li2_tol
  )

  out <- data.frame(
    time = res$grid,
    p_survival = res$p_surv,
    empirical_mean_variance = res$Var_mean_th_cond_surv,
    conditioning = conditioning,
    stringsAsFactors = FALSE
  )
  attr(out, "original_result") <- res
  out
}

simulate_empirical_mean <- function(birth,
                                    death,
                                    sigma2 = 1,
                                    time_start = 0,
                                    time_end,
                                    time_step,
                                    B,
                                    x0 = 0,
                                    seed = NULL,
                                    n_envelope = 4000,
                                    safety_factor = 1.05) {
  requested_grid <- .make_public_grid(time_start, time_end, time_step)
  birth_fun <- .as_rate_function(birth, "birth")
  death_fun <- .as_rate_function(death, "death")
  env <- .bd_env_meanvar()

  raw <- env$simulate_BD_BM_mean_paths(
    lambda_fun = birth_fun,
    mu_fun = death_fun,
    sigma2 = sigma2,
    tmax = time_end,
    dt = time_step,
    B = B,
    x0 = x0,
    seed = seed,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  idx <- .match_grid_indices(raw$grid, requested_grid)

  out <- list(
    time = raw$grid[idx],
    empirical_mean = raw$M[, idx, drop = FALSE],
    n_lineages = raw$N[, idx, drop = FALSE],
    sigma2 = sigma2,
    x0 = x0,
    B = as.integer(B),
    time_start = time_start,
    time_end = time_end,
    time_step = time_step,
    original_result = raw
  )
  out$grid <- out$time
  out$M <- out$empirical_mean
  out$N <- out$n_lineages
  class(out) <- c("empirical_mean_simulation", "list")
  out
}

summarise_empirical_mean_simulation <- function(x) {
  if (!inherits(x, "empirical_mean_simulation")) {
    stop("'x' must be the output of simulate_empirical_mean().", call. = FALSE)
  }
  env <- .bd_env_meanvar()
  raw <- env$compute_empirical_BD_BM_meanvar(.as_source_sim_mean(x))
  data.frame(
    time = raw$grid,
    p_survival_empirical = raw$p_surv_emp,
    empirical_mean_empirical_cond_survival = raw$Mean_emp_cond_surv,
    empirical_mean_variance_empirical_cond_survival = raw$Var_mean_emp_cond_surv
  )
}

plot_empirical_mean_variance <- function(simulation = NULL,
                                         theory = NULL,
                                         summary = NULL,
                                         main = NULL,
                                         xlab = "time",
                                         ylab = "Variance of empirical mean") {
  if (!is.null(simulation) && !inherits(simulation, "empirical_mean_simulation")) {
    stop("'simulation' must be NULL or the output of simulate_empirical_mean().", call. = FALSE)
  }

  if (is.null(summary) && !is.null(simulation)) {
    summary <- summarise_empirical_mean_simulation(simulation)
  }

  y <- c(0)
  if (!is.null(summary)) {
    y <- c(y, summary$empirical_mean_variance_empirical_cond_survival)
  }
  if (!is.null(theory)) {
    y <- c(y, theory$empirical_mean_variance)
  }
  ymax <- max(y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  x_range <- if (!is.null(simulation)) range(simulation$time) else range(theory$time)
  graphics::plot(x_range, c(0, ymax), type = "n", xlab = xlab, ylab = ylab,
                 main = main)

  legend_labels <- character(0)
  legend_lty <- integer(0)
  legend_lwd <- numeric(0)
  legend_col <- character(0)

  if (!is.null(summary)) {
    graphics::lines(summary$time,
                    summary$empirical_mean_variance_empirical_cond_survival,
                    lwd = 3)
    legend_labels <- c(legend_labels, "Monte Carlo")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 3)
    legend_col <- c(legend_col, "black")
  }

  if (!is.null(theory)) {
    graphics::lines(theory$time, theory$empirical_mean_variance,
                    col = "forestgreen", lwd = 3, lty = 2)
    legend_labels <- c(legend_labels, "theory")
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 3)
    legend_col <- c(legend_col, "forestgreen")
  }

  if (length(legend_labels) > 0L) {
    graphics::legend("topleft", legend = legend_labels, lty = legend_lty,
                     lwd = legend_lwd, col = legend_col, bty = "n")
  }

  invisible(list(simulation = simulation, theory = theory, summary = summary))
}

plot_empirical_mean_paths <- function(simulation,
                                      max_paths = 100,
                                      alpha = 0.08,
                                      main = NULL,
                                      xlab = "time",
                                      ylab = "Empirical mean") {
  if (!inherits(simulation, "empirical_mean_simulation")) {
    stop("'simulation' must be the output of simulate_empirical_mean().", call. = FALSE)
  }
  B <- nrow(simulation$empirical_mean)
  path_idx <- seq_len(min(as.integer(max_paths), B))
  yy <- as.numeric(simulation$empirical_mean[path_idx, , drop = FALSE])
  yrange <- range(yy, na.rm = TRUE)
  if (!all(is.finite(yrange)) || diff(yrange) == 0) {
    yrange <- yrange + c(-1, 1)
  }

  graphics::plot(simulation$time, rep(NA_real_, length(simulation$time)),
                 type = "n", xlab = xlab, ylab = ylab, ylim = yrange,
                 main = main)
  col_path <- grDevices::adjustcolor("black", alpha.f = alpha)
  for (b in path_idx) {
    graphics::lines(simulation$time, simulation$empirical_mean[b, ], col = col_path)
  }
  invisible(simulation)
}

run_empirical_mean_variance_experiment <- function(birth,
                                                   death,
                                                   sigma2 = 1,
                                                   time_start = 0,
                                                   time_end,
                                                   time_step,
                                                   B,
                                                   x0 = 0,
                                                   seed = NULL,
                                                   n_steps = 1200,
                                                   n_envelope = 4000,
                                                   safety_factor = 1.05,
                                                   tol_sing = 1e-7,
                                                   li2_tol = 1e-12) {
  sim <- simulate_empirical_mean(
    birth = birth,
    death = death,
    sigma2 = sigma2,
    time_start = time_start,
    time_end = time_end,
    time_step = time_step,
    B = B,
    x0 = x0,
    seed = seed,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  th <- empirical_mean_variance(
    birth = birth,
    death = death,
    sigma2 = sigma2,
    time_start = time_start,
    time_end = time_end,
    time_step = time_step,
    n_steps = n_steps,
    tol_sing = tol_sing,
    li2_tol = li2_tol
  )

  sm <- summarise_empirical_mean_simulation(sim)

  list(simulation = sim, theory = th, summary = sm)
}
