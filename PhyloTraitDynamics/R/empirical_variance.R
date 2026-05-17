## Public API for the expectation of empirical variance under a birth-death
## process with Brownian traits.

empirical_variance_expectation <- function(birth,
                                           death,
                                           sigma2 = 1,
                                           time_start = 0,
                                           time_end,
                                           time_step,
                                           conditioning = c("none", "survival"),
                                           method = c("auto", "numeric", "closed_constant"),
                                           n_steps = 2000,
                                           tol = 1e-8,
                                           constant_tol = 1e-10,
                                           li2_tol = 1e-12) {
  conditioning <- match.arg(conditioning)
  method <- match.arg(method)

  grid <- .make_public_grid(time_start, time_end, time_step)
  birth_fun <- .as_rate_function(birth, "birth")
  death_fun <- .as_rate_function(death, "death")
  env <- .bd_env_empvar()

  res <- env$compute_theory_BD_BM(
    lambda_fun = birth_fun,
    mu_fun = death_fun,
    sigma2 = sigma2,
    grid = grid,
    theory_method = method,
    n_steps = n_steps,
    tol = tol,
    constant_tol = constant_tol,
    li2_tol = li2_tol
  )

  values <- switch(
    conditioning,
    none = res$V_th_uncond,
    survival = res$V_th_cond_same_time
  )

  out <- data.frame(
    time = res$grid,
    p_survival = res$p_surv,
    empirical_variance_expectation = values,
    conditioning = conditioning,
    method = res$theory_method,
    stringsAsFactors = FALSE
  )
  attr(out, "original_result") <- res
  out
}

simulate_empirical_variance <- function(birth,
                                        death,
                                        sigma2 = 1,
                                        time_start = 0,
                                        time_end,
                                        time_step,
                                        B,
                                        seed = NULL,
                                        n_envelope = 4000,
                                        safety_factor = 1.05) {
  requested_grid <- .make_public_grid(time_start, time_end, time_step)
  birth_fun <- .as_rate_function(birth, "birth")
  death_fun <- .as_rate_function(death, "death")
  env <- .bd_env_empvar()

  raw <- env$simulate_BD_BM_paths(
    lambda_fun = birth_fun,
    mu_fun = death_fun,
    sigma2 = sigma2,
    tmax = time_end,
    dt = time_step,
    B = B,
    seed = seed,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  idx <- .match_grid_indices(raw$grid, requested_grid)

  out <- list(
    time = raw$grid[idx],
    empirical_variance = raw$V[, idx, drop = FALSE],
    n_lineages = raw$N[, idx, drop = FALSE],
    sigma2 = sigma2,
    B = as.integer(B),
    time_start = time_start,
    time_end = time_end,
    time_step = time_step,
    original_result = raw
  )
  ## Keep source-compatible names as aliases.  This makes the object usable by
  ## the original summary code without changing the scientific implementation.
  out$grid <- out$time
  out$V <- out$empirical_variance
  out$N <- out$n_lineages
  class(out) <- c("empirical_variance_simulation", "list")
  out
}

summarise_empirical_variance_simulation <- function(x) {
  if (!inherits(x, "empirical_variance_simulation")) {
    stop("'x' must be the output of simulate_empirical_variance().", call. = FALSE)
  }
  env <- .bd_env_empvar()
  raw <- env$compute_empirical_BD_BM(.as_source_sim_empvar(x))
  data.frame(
    time = raw$grid,
    p_survival_empirical = raw$p_surv_emp,
    empirical_variance_expectation_empirical = raw$V_emp_uncond,
    empirical_variance_expectation_empirical_cond_survival = raw$V_emp_cond_same_time
  )
}

plot_empirical_variance_expectation <- function(simulation = NULL,
                                                theory = NULL,
                                                summary = NULL,
                                                show_paths = TRUE,
                                                conditioning = c("none", "survival"),
                                                max_paths = Inf,
                                                alpha = 0.05,
                                                main = NULL,
                                                xlab = "time",
                                                ylab = "Empirical variance") {
  conditioning <- match.arg(conditioning)

  if (!is.null(simulation) && !inherits(simulation, "empirical_variance_simulation")) {
    stop("'simulation' must be NULL or the output of simulate_empirical_variance().", call. = FALSE)
  }

  if (is.null(summary) && !is.null(simulation)) {
    summary <- summarise_empirical_variance_simulation(simulation)
  }

  y <- c(0)
  path_idx <- integer(0)
  if (!is.null(simulation) && show_paths) {
    if (is.finite(max_paths)) {
      path_idx <- seq_len(min(as.integer(max_paths), nrow(simulation$empirical_variance)))
    } else {
      path_idx <- seq_len(nrow(simulation$empirical_variance))
    }
    y <- c(y, as.numeric(simulation$empirical_variance[path_idx, , drop = FALSE]))
  }

  if (!is.null(summary)) {
    if (conditioning == "none") {
      y <- c(y, summary$empirical_variance_expectation_empirical)
    } else {
      y <- c(y, summary$empirical_variance_expectation_empirical_cond_survival)
    }
  }
  if (!is.null(theory)) {
    y <- c(y, theory$empirical_variance_expectation)
  }

  ymax <- max(y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  x_range <- if (!is.null(simulation)) range(simulation$time) else range(theory$time)
  graphics::plot(x_range, c(0, ymax), type = "n", xlab = xlab, ylab = ylab,
                 main = main)

  if (!is.null(simulation) && show_paths && length(path_idx) > 0L) {
    col_path <- grDevices::adjustcolor("black", alpha.f = alpha)
    for (b in path_idx) {
      graphics::lines(simulation$time, simulation$empirical_variance[b, ], col = col_path)
    }
  }

  legend_labels <- character(0)
  legend_lty <- integer(0)
  legend_lwd <- numeric(0)
  legend_col <- character(0)

  if (!is.null(summary)) {
    if (conditioning == "none") {
      yy <- summary$empirical_variance_expectation_empirical
      lab <- "Monte Carlo mean"
    } else {
      yy <- summary$empirical_variance_expectation_empirical_cond_survival
      lab <- "Monte Carlo mean | N(t)>0"
    }
    graphics::lines(summary$time, yy, lwd = 3)
    legend_labels <- c(legend_labels, lab)
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 3)
    legend_col <- c(legend_col, "black")
  }

  if (!is.null(theory)) {
    graphics::lines(theory$time, theory$empirical_variance_expectation,
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

run_empirical_variance_experiment <- function(birth,
                                              death,
                                              sigma2 = 1,
                                              time_start = 0,
                                              time_end,
                                              time_step,
                                              B,
                                              conditioning = c("none", "survival"),
                                              method = c("auto", "numeric", "closed_constant"),
                                              seed = NULL,
                                              n_steps = 2000,
                                              n_envelope = 4000,
                                              safety_factor = 1.05,
                                              tol = 1e-8,
                                              constant_tol = 1e-10,
                                              li2_tol = 1e-12) {
  conditioning <- match.arg(conditioning)
  method <- match.arg(method)

  sim <- simulate_empirical_variance(
    birth = birth,
    death = death,
    sigma2 = sigma2,
    time_start = time_start,
    time_end = time_end,
    time_step = time_step,
    B = B,
    seed = seed,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  th <- empirical_variance_expectation(
    birth = birth,
    death = death,
    sigma2 = sigma2,
    time_start = time_start,
    time_end = time_end,
    time_step = time_step,
    conditioning = conditioning,
    method = method,
    n_steps = n_steps,
    tol = tol,
    constant_tol = constant_tol,
    li2_tol = li2_tol
  )

  sm <- summarise_empirical_variance_simulation(sim)

  list(simulation = sim, theory = th, summary = sm)
}
