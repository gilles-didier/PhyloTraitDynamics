## ================================================================
## Public API for the expectation of empirical variance under a
## birth-death process with Brownian traits.
##
## This file is a package-oriented wrapper around the original script.
## The numerical and simulation routines are kept as internal functions
## with an .empirical_variance_ prefix.
## ================================================================

#' Expected Empirical Variance Under Birth-Death Brownian Dynamics
#'
#' Computes the theoretical expectation of the empirical trait variance on a
#' time grid for a birth-death process with Brownian trait evolution.
#'
#' @param birth Non-negative birth rate, either a numeric constant or a function
#'   of time.
#' @param death Non-negative death rate, either a numeric constant or a function
#'   of time.
#' @param sigma2 Brownian variance parameter.
#' @param time_end,time_step Time grid definition.
#' @param conditioning Conditioning used for the expectation. `"none"` is
#'   unconditional and `"survival"` conditions on `N(t) > 0`.
#' @param method Theory evaluation method. `"auto"` uses the closed constant-rate
#'   formula when the birth and death rates are constant on the time grid, and
#'   otherwise uses the numerical method. `"closed_constant"` requires constant
#'   birth and death rates.
#' @param n_steps Number of integration steps for numerical theory.
#' @param tol,constant_tol,li2_tol Numerical tolerances passed to the original
#'   implementation.
#'
#' @return A data frame with time, survival probability, expected empirical
#'   variance, conditioning, and method columns. The original script result is
#'   attached as attribute `"original_result"`.
#'
#' @export
empirical_variance_compute_expectation <- function(birth,
                                                   death,
                                                   sigma2 = 1,
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

  grid <- .empirical_variance_make_time_grid(time_end, time_step)
  birth_fun <- .empirical_variance_make_rate_function(birth, "birth")
  death_fun <- .empirical_variance_make_rate_function(death, "death")

  res <- .empirical_variance_compute_theory_BD_BM(
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


#' Simulate Empirical Variance Under Birth-Death Brownian Dynamics
#'
#' Simulates birth-death Brownian paths and records the empirical variance of
#' living lineages on the requested time grid.
#'
#' @inheritParams empirical_variance_compute_expectation
#' @param B Number of simulated paths.
#' @param seed Optional random seed.
#' @param n_envelope Number of points used by the thinning envelope in the
#'   original simulator.
#' @param safety_factor Multiplicative safety factor for the thinning envelope.
#'
#' @return An object of class `"empirical_variance_simulation"` containing the
#'   time grid, empirical variances, lineage counts, simulation parameters, and
#'   the original result.
#'
#' @export
empirical_variance_simulate <- function(birth,
                                        death,
                                        sigma2 = 1,
                                        time_end,
                                        time_step,
                                        B,
                                        seed = NULL,
                                        n_envelope = 4000,
                                        safety_factor = 1.05) {
  birth_fun <- .empirical_variance_make_rate_function(birth, "birth")
  death_fun <- .empirical_variance_make_rate_function(death, "death")

  raw <- .empirical_variance_simulate_BD_BM_paths(
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

  out <- list(
    time = raw$grid,
    empirical_variance = raw$V,
    n_lineages = raw$N,
    sigma2 = sigma2,
    B = as.integer(B),
    time_end = time_end,
    time_step = time_step,
    original_result = raw
  )
  out$grid <- out$time
  out$V <- out$empirical_variance
  out$N <- out$n_lineages
  class(out) <- c("empirical_variance_simulation", "list")
  out
}


#' Summarise an Empirical Variance Simulation
#'
#' Computes Monte Carlo summaries from a simulation returned by
#' `empirical_variance_simulate()`.
#'
#' @param x An `"empirical_variance_simulation"` object.
#'
#' @return A data frame with time, empirical survival probability,
#'   unconditional Monte Carlo mean empirical variance, and survival-conditioned
#'   Monte Carlo mean empirical variance.
#'
#' @export
empirical_variance_summarise_simulation <- function(x) {
  if (!inherits(x, "empirical_variance_simulation")) {
    stop("'x' must be the output of empirical_variance_simulate().", call. = FALSE)
  }

  raw <- .empirical_variance_compute_empirical_BD_BM(sim_res = x)

  data.frame(
    time = raw$grid,
    p_survival_empirical = raw$p_surv_emp,
    empirical_variance_expectation_empirical = raw$V_emp_uncond,
    empirical_variance_expectation_empirical_cond_survival = raw$V_emp_cond_same_time
  )
}


#' Plot Expected Empirical Variance
#'
#' Plots simulated empirical variance paths, Monte Carlo summaries, and/or the
#' theoretical expectation on a common set of axes.
#'
#' @param simulation Optional `"empirical_variance_simulation"` object.
#' @param theory Optional data frame returned by
#'   `empirical_variance_compute_expectation()`.
#' @param summary Optional data frame returned by
#'   `empirical_variance_summarise_simulation()`. If omitted and `simulation` is
#'   supplied, it is computed automatically.
#' @param show_paths Logical; whether to draw individual simulated paths.
#' @param conditioning Which expectation to plot.
#' @param max_paths Maximum number of simulated paths to draw.
#' @param alpha Transparency used for simulated paths.
#' @param main,xlab,ylab Base graphics labels.
#'
#' @return Invisibly returns a list containing `simulation`, `theory`, and
#'   `summary`.
#'
#' @export
empirical_variance_plot_expectation <- function(simulation = NULL,
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
    stop("'simulation' must be NULL or the output of empirical_variance_simulate().", call. = FALSE)
  }

  if (is.null(simulation) && is.null(theory)) {
    stop("At least one of 'simulation' and 'theory' must be supplied.", call. = FALSE)
  }

  if (is.null(summary) && !is.null(simulation)) {
    summary <- empirical_variance_summarise_simulation(simulation)
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

  legend_labels <- character(0)
  legend_lty <- integer(0)
  legend_lwd <- numeric(0)
  legend_col <- character(0)

  if (!is.null(simulation) && show_paths && length(path_idx) > 0L) {
    col_path <- grDevices::adjustcolor("black", alpha.f = alpha)
    for (b in path_idx) {
      graphics::lines(simulation$time, simulation$empirical_variance[b, ], col = col_path)
    }
    legend_labels <- c(legend_labels, "simulated trajectories")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 1)
    legend_col <- c(legend_col, col_path)
  }

  if (!is.null(summary)) {
    if (conditioning == "none") {
      yy <- summary$empirical_variance_expectation_empirical
      lab <- "Monte Carlo mean"
    } else {
      yy <- summary$empirical_variance_expectation_empirical_cond_survival
      lab <- "Monte Carlo mean | N(t)>0"
    }
    graphics::lines(summary$time, yy, lwd = 4)
    legend_labels <- c(legend_labels, lab)
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 4)
    legend_col <- c(legend_col, "black")
  }

  if (!is.null(theory)) {
    graphics::lines(theory$time, theory$empirical_variance_expectation,
                    col = "green", lwd = 3, lty = 2)
    legend_labels <- c(legend_labels, "expectation")
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 3)
    legend_col <- c(legend_col, "green")
  }

  if (length(legend_labels) > 0L) {
    graphics::legend("topleft", legend = legend_labels, lty = legend_lty,
                     lwd = legend_lwd, col = legend_col, bty = "n")
  }

  invisible(list(simulation = simulation, theory = theory, summary = summary))
}


.empirical_variance_make_rate_function <- function(rate, name = "rate") {
  if (is.function(rate)) {
    return(function(t) {
      y <- try(rate(t), silent = TRUE)

      ## If the function is not vectorized, evaluate pointwise.
      if (inherits(y, "try-error") || length(y) != length(t)) {
        y <- vapply(t, rate, numeric(1))
      }

      y <- as.numeric(y)

      if (length(y) != length(t)) {
        stop(name, " must return either one value per input time or be scalar.")
      }

      y
    })
  }

  if (is.numeric(rate) && length(rate) == 1L) {
    return(function(t) rep(rate, length(t)))
  }

  stop(name, " must be either a function of time or a single numeric value.")
}


## ================================================================
## 0) Generic utilities
## ================================================================

.empirical_variance_make_time_grid <- function(tmax, dt){
  if (!is.finite(tmax) || tmax < 0){
    stop("'tmax' must be a finite nonnegative number.")
  }
  if (!is.finite(dt) || dt <= 0){
    stop("'dt' must be a strictly positive number.")
  }

  grid <- seq(0, tmax, by = dt)
  if (length(grid) == 0L || abs(tail(grid, 1L) - tmax) > 1e-14){
    grid <- c(grid, tmax)
  }
  unique(grid)
}

.empirical_variance_trapz_vec <- function(x, y){
  n <- length(x)
  if (n != length(y)){
    stop("'x' and 'y' must have the same length.")
  }
  if (n <= 1L){
    return(0)
  }
  sum(diff(x) * (y[-n] + y[-1L]) / 2)
}

.empirical_variance_cumtrapz_vec <- function(x, y){
  n <- length(x)
  if (n != length(y)){
    stop("'x' and 'y' must have the same length.")
  }
  out <- numeric(n)
  if (n <= 1L){
    return(out)
  }
  out[-1L] <- cumsum(diff(x) * (y[-n] + y[-1L]) / 2)
  out
}

.empirical_variance_eval_on_grid_scalar <- function(x, y, xout){
  approx(x = x, y = y, xout = xout, rule = 2)$y
}

.empirical_variance_eval_rate_fun <- function(rate_fun, x){
  out <- rate_fun(x)
  if (length(out) == 1L && length(x) > 1L){
    out <- rep(out, length(x))
  }
  if (length(out) != length(x)){
    stop("A rate function must return either a scalar or a vector of the same length as its input.")
  }
  as.numeric(out)
}

.empirical_variance_check_rate_functions <- function(lambda_fun, mu_fun, grid){
  lambda_vals <- .empirical_variance_eval_rate_fun(lambda_fun, grid)
  mu_vals <- .empirical_variance_eval_rate_fun(mu_fun, grid)

  if (any(!is.finite(lambda_vals)) || any(!is.finite(mu_vals))){
    stop("Non-finite rate values detected.")
  }
  if (any(lambda_vals < 0) || any(mu_vals < 0)){
    stop("Rates must be nonnegative on the working grid.")
  }

  invisible(list(lambda = lambda_vals, mu = mu_vals))
}

.empirical_variance_poly_rate_fun <- function(coeffs, floor_value = 0){
  force(coeffs)
  force(floor_value)

  function(t){
    t <- as.numeric(t)
    out <- rep(0, length(t))
    for (a in rev(coeffs)){
      out <- out * t + a
    }
    pmax(floor_value, out)
  }
}


.empirical_variance_constant_rate_fun <- function(rate){
  if (!is.finite(rate) || rate < 0){
    stop("'rate' must be a finite nonnegative number.")
  }
  force(rate)
  function(t){
    rep(rate, length(as.numeric(t)))
  }
}

.empirical_variance_is_constant_on_grid <- function(values, tol = 1e-10){
  values <- as.numeric(values)
  if (length(values) == 0L || any(!is.finite(values))){
    return(FALSE)
  }
  max(abs(values - values[1L])) <= tol * max(1, abs(values[1L]))
}

.empirical_variance_detect_constant_rates <- function(lambda_fun, mu_fun, grid, tol = 1e-10){
  lambda_vals <- .empirical_variance_eval_rate_fun(lambda_fun, grid)
  mu_vals <- .empirical_variance_eval_rate_fun(mu_fun, grid)

  if (!.empirical_variance_is_constant_on_grid(lambda_vals, tol = tol) ||
      !.empirical_variance_is_constant_on_grid(mu_vals, tol = tol)){
    return(NULL)
  }

  lambda <- lambda_vals[1L]
  mu <- mu_vals[1L]

  if (lambda < 0 || mu < 0){
    stop("Rates must be nonnegative.")
  }

  list(lambda = lambda, mu = mu)
}

.empirical_variance_li2_series_0_05 <- function(x, tol = 1e-12, max_iter = 10000L){
  ## Power series for |x| <= 1/2:
  ## Li_2(x) = sum_{k >= 1} x^k / k^2.
  if (x == 0){
    return(0)
  }

  term <- x
  out <- term

  for (k in 2:max_iter){
    term <- term * x
    add <- term / (k * k)
    out <- out + add

    if (abs(add) <= tol * max(1, abs(out))){
      return(out)
    }
  }

  warning("Dilogarithm series did not reach the requested tolerance.")
  out
}

.empirical_variance_li2_real_scalar <- function(z, tol = 1e-12, max_iter = 10000L){
  ## Real dilogarithm Li_2(z) on the real domain z <= 1.
  ## This is enough for the constant birth-death closed forms below.
  if (!is.finite(z)){
    stop("Non-finite argument passed to .empirical_variance_li2_real().")
  }
  if (z > 1 + 1e-12){
    stop(".empirical_variance_li2_real() is implemented only for real arguments <= 1.")
  }
  if (z > 1){
    z <- 1
  }
  if (z == 1){
    return(pi^2 / 6)
  }
  if (z == 0){
    return(0)
  }

  if (z < 0){
    ## Li_2(z) = -Li_2(z/(z-1)) - 1/2 log(1-z)^2.
    w <- z / (z - 1)
    return(-.empirical_variance_li2_real_scalar(w, tol = tol, max_iter = max_iter) -
             0.5 * log1p(-z)^2)
  }

  if (z <= 0.5){
    return(.empirical_variance_li2_series_0_05(z, tol = tol, max_iter = max_iter))
  }

  ## Li_2(z) = pi^2/6 - log(z) log(1-z) - Li_2(1-z).
  pi^2 / 6 - log(z) * log1p(-z) -
    .empirical_variance_li2_real_scalar(1 - z, tol = tol, max_iter = max_iter)
}

.empirical_variance_li2_real <- function(z, tol = 1e-12, max_iter = 10000L){
  vapply(
    as.numeric(z),
    .empirical_variance_li2_real_scalar,
    numeric(1),
    tol = tol,
    max_iter = as.integer(max_iter)
  )
}

.empirical_variance_constant_bd_survival <- function(lambda, mu, t, critical_tol = 1e-10){
  t <- as.numeric(t)
  if (lambda < 0 || mu < 0){
    stop("Rates must be nonnegative.")
  }

  if (abs(lambda - mu) <= critical_tol * max(1, abs(lambda), abs(mu))){
    if (lambda == 0){
      return(rep(1, length(t)))
    }
    return(1 / (1 + lambda * t))
  }

  r <- lambda - mu
  r / (lambda - mu * exp(-r * t))
}

.empirical_variance_constant_bd_phi <- function(lambda, mu, t, T = t, critical_tol = 1e-10){
  ## Time change phi_T(t) for constant rates.
  t <- as.numeric(t)
  if (lambda < 0 || mu < 0){
    stop("Rates must be nonnegative.")
  }

  if (abs(lambda - mu) <= critical_tol * max(1, abs(lambda), abs(mu))){
    if (lambda == 0){
      return(rep(0, length(t)))
    }
    return(log((1 + lambda * T) / (1 + lambda * (T - t))))
  }

  r <- lambda - mu
  r * t - log((lambda - mu * exp(-r * (T - t))) /
                (lambda - mu * exp(-r * T)))
}

.empirical_variance_closed_constant_empvar_uncond <- function(lambda, mu, sigma2, t,
                                          critical_tol = 1e-10,
                                          li2_tol = 1e-12){
  t <- as.numeric(t)
  out <- numeric(length(t))

  if (lambda < 0 || mu < 0){
    stop("Rates must be nonnegative.")
  }
  if (!is.finite(sigma2) || sigma2 < 0){
    stop("'sigma2' must be a finite nonnegative number.")
  }
  if (sigma2 == 0 || lambda == 0){
    return(out)
  }

  is_zero_t <- abs(t) <= 1e-14
  if (all(is_zero_t)){
    return(out)
  }

  nz <- !is_zero_t
  tt <- t[nz]

  if (abs(lambda - mu) <= critical_tol * max(1, abs(lambda), abs(mu))){
    z <- lambda * tt / (1 + lambda * tt)
    bracket <-
      lambda * tt^2 / (1 + lambda * tt) -
      4 * tt +
      (2 / lambda) * log1p(lambda * tt) +
      (2 * (1 + lambda * tt) / lambda) * .empirical_variance_li2_real(z, tol = li2_tol)

    out[nz] <- sigma2 / (1 + lambda * tt) * bracket
    return(pmax(0, out))
  }

  r <- lambda - mu
  e <- exp(-r * tt)

  if (mu == 0){
    ## Yule limit of the non-critical formula. Note the + e^{-lambda t} Li_2 term
    ## inside the parenthesis below; this is the sign obtained by setting mu=0
    ## in the general non-critical expression.
    e_yule <- exp(-lambda * tt)
    z <- 1 - exp(lambda * tt)
    out[nz] <- sigma2 * (
      tt * (1 - e_yule) -
        (2 / lambda) * (1 - e_yule + e_yule * .empirical_variance_li2_real(z, tol = li2_tol))
    )
    return(pmax(0, out))
  }

  denom <- lambda - mu * e
  log_arg <- denom / (lambda - mu)
  li2_arg_1 <- 1 - exp(r * tt)
  li2_arg_2 <- mu * (1 - e) / denom

  if (any(log_arg <= 0) || any(li2_arg_1 > 1 + 1e-10) || any(li2_arg_2 > 1 + 1e-10)){
    stop("Closed-form expression received arguments outside their real domain.")
  }

  bracket <-
    tt * (
      lambda * (1 - e) / denom -
        (2 * mu / lambda) * e
    ) -
    2 * (1 - e) / r +
    (2 * e / lambda) * log(log_arg) -
    (2 * e * denom / (lambda * r)) *
    (.empirical_variance_li2_real(li2_arg_1, tol = li2_tol) -
       .empirical_variance_li2_real(li2_arg_2, tol = li2_tol))

  out[nz] <- sigma2 * r / denom * bracket
  pmax(0, out)
}

.empirical_variance_compute_theory_BD_BM_closed_constant <- function(lambda, mu, sigma2, grid,
                                                 critical_tol = 1e-10,
                                                 li2_tol = 1e-12){
  if (any(!is.finite(grid)) || any(grid < 0)){
    stop("'grid' must contain finite nonnegative times.")
  }
  if (is.unsorted(grid, strictly = FALSE)){
    stop("'grid' must be sorted increasingly.")
  }

  p_surv <- .empirical_variance_constant_bd_survival(
    lambda = lambda,
    mu = mu,
    t = grid,
    critical_tol = critical_tol
  )

  V_uncond <- .empirical_variance_closed_constant_empvar_uncond(
    lambda = lambda,
    mu = mu,
    sigma2 = sigma2,
    t = grid,
    critical_tol = critical_tol,
    li2_tol = li2_tol
  )

  V_cond_same_time <- rep(NA_real_, length(grid))
  keep <- p_surv > 0
  V_cond_same_time[keep] <- V_uncond[keep] / p_surv[keep]

  list(
    grid = grid,
    sigma2 = sigma2,
    theory_method = "closed_constant",
    lambda_constant = lambda,
    mu_constant = mu,
    p_surv = p_surv,
    V_th_uncond = V_uncond,
    V_th_cond_same_time = V_cond_same_time,
    reconstructed = NULL
  )
}


## ================================================================
## 1) Reconstructed process at fixed final time t
##    Solve backward:
##      g_t'(w) = (lambda(w) - mu(w)) g_t(w) - lambda(w),
##      g_t(t)  = 1.
## ================================================================

.empirical_variance_solve_reconstructed_process <- function(lambda_fun, mu_fun, t,
                                        n_steps = 2000,
                                        check_rates = TRUE){
  if (!is.finite(t) || t < 0){
    stop("'t' must be a finite nonnegative number.")
  }
  if (!is.numeric(n_steps) || length(n_steps) != 1L || n_steps < 1){
    stop("'n_steps' must be a positive integer.")
  }
  n_steps <- as.integer(n_steps)

  if (t == 0){
    grid <- 0
    lambda0 <- .empirical_variance_eval_rate_fun(lambda_fun, 0)
    mu0 <- .empirical_variance_eval_rate_fun(mu_fun, 0)

    if (!is.finite(lambda0) || !is.finite(mu0) || lambda0 < 0 || mu0 < 0){
      stop("Invalid rate values at time 0.")
    }

    return(list(
      t = t,
      grid = grid,
      lambda = lambda0,
      mu = mu0,
      g = 1,
      p_surv = 1,
      p_extinct = 0,
      rebirth = lambda0,
      phi = 0
    ))
  }

  grid_desc <- seq(t, 0, length.out = n_steps + 1L)

  if (check_rates){
    .empirical_variance_check_rate_functions(lambda_fun, mu_fun, grid_desc)
  }

  rhs_g <- function(w, g){
    (.empirical_variance_eval_rate_fun(lambda_fun, w) - .empirical_variance_eval_rate_fun(mu_fun, w)) * g -
      .empirical_variance_eval_rate_fun(lambda_fun, w)
  }

  g_desc <- numeric(length(grid_desc))
  g_desc[1L] <- 1

  for (i in seq_len(n_steps)){
    w <- grid_desc[i]
    h <- grid_desc[i + 1L] - grid_desc[i]
    g_i <- g_desc[i]

    k1 <- rhs_g(w, g_i)
    k2 <- rhs_g(w + h / 2, g_i + h * k1 / 2)
    k3 <- rhs_g(w + h / 2, g_i + h * k2 / 2)
    k4 <- rhs_g(w + h, g_i + h * k3)

    g_next <- g_i + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    if (!is.finite(g_next) || g_next <= 0){
      stop(
        paste0(
          "The numerical solution for g became nonpositive or non-finite. ",
          "Try increasing 'n_steps' or using smoother/smaller rates."
        )
      )
    }

    g_desc[i + 1L] <- g_next
  }

  grid <- rev(grid_desc)
  g <- rev(g_desc)
  lambda_vals <- .empirical_variance_eval_rate_fun(lambda_fun, grid)
  mu_vals <- .empirical_variance_eval_rate_fun(mu_fun, grid)

  rebirth <- lambda_vals / g
  phi <- .empirical_variance_cumtrapz_vec(grid, rebirth)

  list(
    t = t,
    grid = grid,
    lambda = lambda_vals,
    mu = mu_vals,
    g = g,
    p_surv = 1 / g,
    p_extinct = 1 - 1 / g,
    rebirth = rebirth,
    phi = phi
  )
}


## ================================================================
## 2) Stable integrands for the theoretical expectations
## ================================================================

.empirical_variance_fprob_from_phi <- function(phi_s, phi_t, tol = 1e-8){
  e_minus_phi_t <- exp(-phi_t)
  out <- numeric(length(phi_s))

  small <- abs(phi_s) < tol
  large <- !small

  if (any(large)){
    x <- phi_s[large]
    exm1 <- expm1(x)
    ex <- 1 + exm1
    num <- (1 - e_minus_phi_t * ex) * (1 - ex + x * ex)
    out[large] <- 1 - e_minus_phi_t - 2 * num / (exm1^2)
  }

  if (any(small)){
    x <- phi_s[small]
    out[small] <- ((1 + 2 * e_minus_phi_t) / 3) * x +
      (e_minus_phi_t / 6) * x^2 +
      ((e_minus_phi_t - 1) / 90) * x^3
  }

  out
}

## ================================================================
## 3) Theory on a time grid
## ================================================================

.empirical_variance_compute_theory_BD_BM <- function(lambda_fun, mu_fun, sigma2, grid,
                                 theory_method = c("numeric", "closed_constant", "auto"),
                                 n_steps = 2000,
                                 tol = 1e-8,
                                 constant_tol = 1e-10,
                                 li2_tol = 1e-12,
                                 return_reconstructed = FALSE){
  if (any(!is.finite(grid)) || any(grid < 0)){
    stop("'grid' must contain finite nonnegative times.")
  }
  if (is.unsorted(grid, strictly = FALSE)){
    stop("'grid' must be sorted increasingly.")
  }
  if (!is.finite(sigma2) || sigma2 < 0){
    stop("'sigma2' must be a finite nonnegative number.")
  }

  theory_method <- match.arg(theory_method)

  constant_rates <- .empirical_variance_detect_constant_rates(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    grid = grid,
    tol = constant_tol
  )

  if (theory_method == "auto"){
    theory_method <- if (is.null(constant_rates)) "numeric" else "closed_constant"
  }

  if (theory_method == "closed_constant"){
    if (is.null(constant_rates)){
      stop("'theory_method = closed_constant' requires rates that are constant on the time grid.")
    }
    return(.empirical_variance_compute_theory_BD_BM_closed_constant(
      lambda = constant_rates$lambda,
      mu = constant_rates$mu,
      sigma2 = sigma2,
      grid = grid,
      critical_tol = constant_tol,
      li2_tol = li2_tol
    ))
  }

  n_grid <- length(grid)
  V_uncond <- numeric(n_grid)
  V_cond_same_time <- numeric(n_grid)
  p_surv <- numeric(n_grid)

  reconstructed <- if (return_reconstructed) vector("list", n_grid) else NULL

  for (i in seq_along(grid)){
    t <- grid[i]

    rec_t <- .empirical_variance_solve_reconstructed_process(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      t = t,
      n_steps = n_steps
    )

    phi_t <- tail(rec_t$phi, 1L)
    f_vals <- .empirical_variance_fprob_from_phi(rec_t$phi, phi_t, tol = tol)
    int_f <- .empirical_variance_trapz_vec(rec_t$grid, f_vals)

    p_surv[i] <- rec_t$p_surv[1L]
    V_cond_same_time[i] <- sigma2 * int_f
    V_uncond[i] <- sigma2 * p_surv[i] * int_f

    if (return_reconstructed){
      reconstructed[[i]] <- rec_t
    }
  }

  list(
    grid = grid,
    sigma2 = sigma2,
    theory_method = "numeric",
    p_surv = p_surv,
    V_th_uncond = V_uncond,
    V_th_cond_same_time = V_cond_same_time,
    reconstructed = reconstructed
  )
}

## ================================================================
## 4) Inhomogeneous birth-death simulation by thinning
## ================================================================

.empirical_variance_make_rate_envelope <- function(lambda_fun, mu_fun, tmax,
                               n_envelope = 4000,
                               safety_factor = 1.05,
                               use_midpoints = TRUE){
  if (!is.numeric(n_envelope) || length(n_envelope) != 1L || n_envelope < 10){
    stop("'n_envelope' must be an integer >= 10.")
  }
  if (!is.finite(safety_factor) || safety_factor < 1){
    stop("'safety_factor' must be a finite number >= 1.")
  }
  n_envelope <- as.integer(n_envelope)

  coarse_grid <- seq(0, tmax, length.out = n_envelope + 1L)
  if (isTRUE(use_midpoints) && length(coarse_grid) >= 2L){
    mid_grid <- (coarse_grid[-1L] + coarse_grid[-length(coarse_grid)]) / 2
    grid <- sort(unique(c(coarse_grid, mid_grid)))
  } else {
    grid <- coarse_grid
  }

  checked <- .empirical_variance_check_rate_functions(lambda_fun, mu_fun, grid)
  rate_sum <- checked$lambda + checked$mu
  suffix_max <- safety_factor * rev(cummax(rev(rate_sum)))

  upper_from <- function(time){
    idx <- findInterval(time, grid, rightmost.closed = TRUE, all.inside = TRUE)
    suffix_max[idx]
  }

  list(
    grid = grid,
    rate_sum = rate_sum,
    suffix_max = suffix_max,
    safety_factor = safety_factor,
    upper_from = upper_from
  )
}

.empirical_variance_sample_next_event_thinning <- function(time, k, tmax,
                                       lambda_fun, mu_fun,
                                       envelope,
                                       max_restart = 20L){
  if (k <= 0 || time >= tmax){
    return(Inf)
  }
  if (!is.numeric(max_restart) || length(max_restart) != 1L || max_restart < 1){
    stop("'max_restart' must be a positive integer.")
  }
  max_restart <- as.integer(max_restart)

  local_upper <- envelope$upper_from(time)
  if (!is.finite(local_upper) || local_upper < 0){
    stop("Invalid upper bound in thinning.")
  }
  if (local_upper == 0){
    return(Inf)
  }

  safety_factor <- if (!is.null(envelope$safety_factor)) envelope$safety_factor else 1.05

  for (restart in seq_len(max_restart)){
    M <- k * local_upper
    current <- time

    repeat {
      current <- current + rexp(1, rate = M)
      if (current > tmax){
        return(Inf)
      }

      lambda_now <- .empirical_variance_eval_rate_fun(lambda_fun, current)
      mu_now <- .empirical_variance_eval_rate_fun(mu_fun, current)
      rate_now_one_lineage <- lambda_now + mu_now
      rate_now <- k * rate_now_one_lineage

      if (!is.finite(rate_now) || rate_now < 0){
        stop("Invalid instantaneous rate during thinning.")
      }

      if (rate_now > M * (1 + 1e-12)){
        local_upper <- max(local_upper, as.numeric(rate_now_one_lineage) * safety_factor)
        break
      }

      if (runif(1) <= rate_now / M){
        return(current)
      }
    }
  }

  stop(
    paste0(
      "The thinning envelope remained too small after ", max_restart,
      " adaptive restarts. Increase 'n_envelope' or 'safety_factor'."
    )
  )
}

.empirical_variance_sim_BD_BM_var_path_inhom <- function(lambda_fun, mu_fun, sigma2,
                                     tmax, dt,
                                     envelope = NULL,
                                     n_envelope = 4000,
                                     safety_factor = 1.05){
  if (!is.finite(sigma2) || sigma2 < 0){
    stop("'sigma2' must be a finite nonnegative number.")
  }

  grid <- .empirical_variance_make_time_grid(tmax, dt)
  n_grid <- length(grid)

  if (is.null(envelope)){
    envelope <- .empirical_variance_make_rate_envelope(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      tmax = tmax,
      n_envelope = n_envelope,
      safety_factor = safety_factor
    )
  }

  time <- 0
  traits <- 0
  k <- 1L
  idx <- 1L

  var_path <- numeric(n_grid)
  n_path <- integer(n_grid)
  var_path[1L] <- 0
  n_path[1L] <- 1L

  next_event <- .empirical_variance_sample_next_event_thinning(
    time = time,
    k = k,
    tmax = tmax,
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    envelope = envelope
  )

  while (idx < n_grid){
    next_grid_time <- grid[idx + 1L]
    t_next <- min(next_event, next_grid_time, tmax)

    dt_seg <- t_next - time
    if (dt_seg > 0 && k > 0){
      traits <- traits + rnorm(k, mean = 0, sd = sqrt(sigma2 * dt_seg))
    }
    time <- t_next

    if (abs(time - next_grid_time) < 1e-12){
      idx <- idx + 1L
      n_path[idx] <- k
      var_path[idx] <- if (k >= 2L) var(traits) else 0
      next
    }

    if (abs(time - next_event) < 1e-12){
      lambda_now <- .empirical_variance_eval_rate_fun(lambda_fun, time)
      mu_now <- .empirical_variance_eval_rate_fun(mu_fun, time)
      total_now <- lambda_now + mu_now

      if (total_now <= 0){
        next_event <- Inf
        next
      }

      if (runif(1) < lambda_now / total_now){
        parent <- sample.int(k, 1L)
        traits <- c(traits, traits[parent])
        k <- k + 1L
      } else {
        victim <- sample.int(k, 1L)
        if (k == 1L){
          traits <- numeric(0)
          k <- 0L
        } else {
          traits <- traits[-victim]
          k <- k - 1L
        }
      }

      next_event <- .empirical_variance_sample_next_event_thinning(
        time = time,
        k = k,
        tmax = tmax,
        lambda_fun = lambda_fun,
        mu_fun = mu_fun,
        envelope = envelope
      )

      next
    }

    break
  }

  list(
    grid = grid,
    var = var_path,
    n = n_path
  )
}

.empirical_variance_simulate_BD_BM_paths <- function(lambda_fun, mu_fun, sigma2,
                                 tmax, dt, B,
                                 seed = NULL,
                                 n_envelope = 4000,
                                 safety_factor = 1.05){
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 1){
    stop("'B' must be a positive integer.")
  }
  B <- as.integer(B)

  grid <- .empirical_variance_make_time_grid(tmax, dt)
  n_grid <- length(grid)

  envelope <- .empirical_variance_make_rate_envelope(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    tmax = tmax,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  V <- matrix(0, nrow = B, ncol = n_grid)
  N <- matrix(0L, nrow = B, ncol = n_grid)

  for (b in seq_len(B)){
    out <- .empirical_variance_sim_BD_BM_var_path_inhom(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      sigma2 = sigma2,
      tmax = tmax,
      dt = dt,
      envelope = envelope,
      safety_factor = safety_factor
    )
    V[b, ] <- out$var
    N[b, ] <- out$n
  }

  list(
    grid = grid,
    sigma2 = sigma2,
    V = V,
    N = N,
    B = B,
    tmax = tmax,
    dt = dt
  )
}


## ================================================================
## 5) Empirical summaries extracted from simulations
## ================================================================

.empirical_variance_compute_empirical_BD_BM <- function(sim_res){
  grid <- sim_res$grid
  V <- sim_res$V
  N <- sim_res$N

  ## V is already 0 when N(t) < 2. Thus V_emp_uncond estimates
  ## E[S^2(t)] with this convention.
  alive_same_time <- N > 0
  p_surv_emp <- colMeans(alive_same_time)
  V_emp_uncond <- colMeans(V)

  ## This estimates E[S^2(t) | N(t) > 0], still with S^2(t)=0
  ## when N(t)=1. There is no conditioning on N(t) >= 2.
  V_emp_cond_same_time <- rep(NA_real_, length(grid))
  for (j in seq_along(grid)){
    keep <- alive_same_time[, j]
    if (any(keep)){
      V_emp_cond_same_time[j] <- mean(V[keep, j])
    }
  }

  list(
    grid = grid,
    p_surv_emp = p_surv_emp,
    V_emp_uncond = V_emp_uncond,
    V_emp_cond_same_time = V_emp_cond_same_time
  )
}

## ================================================================
## 6) Comparison helpers
## ================================================================

.empirical_variance_compare_BD_BM_fit <- function(sim_res, theory_res = NULL, empirical_res = NULL){
  if (is.null(empirical_res)){
    empirical_res <- .empirical_variance_compute_empirical_BD_BM(sim_res = sim_res)
  }

  out <- data.frame(
    curve = c("unconditional", "conditional_same_time", "survival_probability"),
    max_abs_diff = NA_real_,
    rmse = NA_real_,
    final_emp = c(
      tail(empirical_res$V_emp_uncond, 1L),
      tail(empirical_res$V_emp_cond_same_time, 1L),
      tail(empirical_res$p_surv_emp, 1L)
    ),
    final_th = NA_real_
  )

  if (is.null(theory_res)){
    return(out)
  }

  diff_uncond <- empirical_res$V_emp_uncond - theory_res$V_th_uncond
  diff_cond_same <- empirical_res$V_emp_cond_same_time - theory_res$V_th_cond_same_time
  diff_surv <- empirical_res$p_surv_emp - theory_res$p_surv

  out$max_abs_diff[1L] <- max(abs(diff_uncond), na.rm = TRUE)
  out$rmse[1L] <- sqrt(mean(diff_uncond^2, na.rm = TRUE))
  out$final_th[1L] <- tail(theory_res$V_th_uncond, 1L)

  out$max_abs_diff[2L] <- max(abs(diff_cond_same), na.rm = TRUE)
  out$rmse[2L] <- sqrt(mean(diff_cond_same^2, na.rm = TRUE))
  out$final_th[2L] <- tail(theory_res$V_th_cond_same_time, 1L)

  out$max_abs_diff[3L] <- max(abs(diff_surv), na.rm = TRUE)
  out$rmse[3L] <- sqrt(mean(diff_surv^2, na.rm = TRUE))
  out$final_th[3L] <- tail(theory_res$p_surv, 1L)

  out
}

## ================================================================
## 7) Plotting
##    This function takes simulation output (and optionally theory)
## ================================================================

.empirical_variance_plot_BD_BM_results <- function(sim_res,
                               theory_res = NULL,
                               empirical_res = NULL,
                               show_paths = TRUE,
                               show_unconditional = TRUE,
                               show_cond_same_time = FALSE,
                               alpha = 0.05,
                               lwd_paths = 1,
                               max_paths = Inf,
                               main = NULL,
                               xlab = "time",
                               ylab = "Empirical variance"){
  grid <- sim_res$grid
  V <- sim_res$V

  if (is.null(empirical_res)){
    empirical_res <- .empirical_variance_compute_empirical_BD_BM(sim_res = sim_res)
  }

  y_candidates <- c(0)
  if (show_paths){
    if (is.finite(max_paths)){
      path_idx <- seq_len(min(as.integer(max_paths), nrow(V)))
    } else {
      path_idx <- seq_len(nrow(V))
    }
    y_candidates <- c(y_candidates, V[path_idx, , drop = FALSE])
  } else {
    path_idx <- integer(0)
  }
  if (show_unconditional){
    y_candidates <- c(y_candidates, empirical_res$V_emp_uncond)
    if (!is.null(theory_res)){
      y_candidates <- c(y_candidates, theory_res$V_th_uncond)
    }
  }
  if (show_cond_same_time){
    y_candidates <- c(y_candidates, empirical_res$V_emp_cond_same_time)
    if (!is.null(theory_res)){
      y_candidates <- c(y_candidates, theory_res$V_th_cond_same_time)
    }
  }

  ymax <- max(y_candidates, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0){
    ymax <- 1
  }

  plot(
    grid, rep(NA_real_, length(grid)),
    type = "n",
    xlab = xlab,
    ylab = ylab,
    ylim = c(0, ymax),
    main = main
  )

  col_path <- adjustcolor("black", alpha.f = alpha)
  if (show_paths && length(path_idx) > 0L){
    for (b in path_idx){
      lines(grid, V[b, ], col = col_path, lwd = lwd_paths)
    }
  }

  legend_labels <- character(0)
  legend_lty <- integer(0)
  legend_lwd <- numeric(0)
  legend_col <- character(0)

  if (show_paths){
    legend_labels <- c(legend_labels, "simulated trajectories")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, lwd_paths)
    legend_col <- c(legend_col, col_path)
  }

  if (show_unconditional){
    lines(grid, empirical_res$V_emp_uncond, lwd = 4)
    legend_labels <- c(legend_labels, "simulations mean")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 4)
    legend_col <- c(legend_col, "black")

    if (!is.null(theory_res)){
#      theory_label <- if (!is.null(theory_res$theory_method)) {
#        paste0("theory, ", theory_res$theory_method)
#      } else {
#        "theory"
#      }
      theory_label <- "expectation"
      lines(grid, theory_res$V_th_uncond, col = "green", lwd = 3, lty = 2)
      legend_labels <- c(legend_labels, theory_label)
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 3)
      legend_col <- c(legend_col, "green")
    }
  }

  if (show_cond_same_time){
    lines(grid, empirical_res$V_emp_cond_same_time, col = "black", lwd = 4)
    legend_labels <- c(legend_labels, "simulations mean")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 4)
    legend_col <- c(legend_col, "black")

    if (!is.null(theory_res)){
      lines(grid, theory_res$V_th_cond_same_time, col = "green", lwd = 3, lty = 2)
      legend_labels <- c(legend_labels, "expectation")
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 3)
      legend_col <- c(legend_col, "green")
    }
  }

  legend(
    "topleft",
    legend = legend_labels,
    lty = legend_lty,
    lwd = legend_lwd,
    col = legend_col,
    bty = "n"
  )

  invisible(list(sim = sim_res, empirical = empirical_res, theory = theory_res))
}

