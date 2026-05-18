## ================================================================
## Uses shared birth-death utilities from birth_death_utils.R.
## ================================================================

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

  grid <- .birth_death_make_time_grid(time_end, time_step)
  birth_fun <- .birth_death_make_rate_function(birth, "birth")
  death_fun <- .birth_death_make_rate_function(death, "death")

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
  birth_fun <- .birth_death_make_rate_function(birth, "birth")
  death_fun <- .birth_death_make_rate_function(death, "death")

  raw <- .birth_death_brownian_simulate_paths(
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
    empirical_mean = raw$M,
    empirical_variance = raw$V,
    n_lineages = raw$N,
    sigma2 = sigma2,
    x0 = 0,
    B = as.integer(B),
    time_end = time_end,
    time_step = time_step,
    original_result = raw
  )
  out$grid <- out$time
  out$M <- out$empirical_mean
  out$V <- out$empirical_variance
  out$N <- out$n_lineages
  out$summary <- .birth_death_brownian_summarise_replicates(out)
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
#' @keywords internal
empirical_variance_summarise_simulation <- function(x) {
  if (!inherits(x, "empirical_variance_simulation")) {
    stop("'x' must be the output of empirical_variance_simulate().", call. = FALSE)
  }

  if (!is.null(x$summary)) {
    return(x$summary)
  }

  .birth_death_brownian_summarise_replicates(x)
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


#' Plot Empirical Variance Paths
#'
#' Plots individual simulated empirical variance paths.
#'
#' @param simulation An `"empirical_variance_simulation"` object.
#' @param max_paths Maximum number of simulated paths to draw.
#' @param alpha Transparency used for simulated paths.
#' @param main,xlab,ylab Base graphics labels.
#'
#' @return Invisibly returns `simulation`.
#'
#' @export
empirical_variance_plot_paths <- function(simulation,
                                          max_paths = 100,
                                          alpha = 0.08,
                                          main = NULL,
                                          xlab = "time",
                                          ylab = "Empirical variance") {
  if (!inherits(simulation, "empirical_variance_simulation")) {
    stop("'simulation' must be the output of empirical_variance_simulate().", call. = FALSE)
  }

  B <- nrow(simulation$empirical_variance)
  path_idx <- seq_len(min(as.integer(max_paths), B))
  yy <- as.numeric(simulation$empirical_variance[path_idx, , drop = FALSE])
  yrange <- range(yy, na.rm = TRUE)
  if (!all(is.finite(yrange)) || diff(yrange) == 0) {
    yrange <- yrange + c(-1, 1)
  }

  graphics::plot(simulation$time, rep(NA_real_, length(simulation$time)),
                 type = "n", xlab = xlab, ylab = ylab, ylim = yrange,
                 main = main)
  col_path <- grDevices::adjustcolor("black", alpha.f = alpha)
  for (b in path_idx) {
    graphics::lines(simulation$time, simulation$empirical_variance[b, ], col = col_path)
  }

  invisible(simulation)
}


## ================================================================
## 0) Generic utilities
## ================================================================

.empirical_variance_eval_on_grid_scalar <- function(x, y, xout){
  approx(x = x, y = y, xout = xout, rule = 2)$y
}

.empirical_variance_is_constant_on_grid <- function(values, tol = 1e-10){
  values <- as.numeric(values)
  if (length(values) == 0L || any(!is.finite(values))){
    return(FALSE)
  }
  max(abs(values - values[1L])) <= tol * max(1, abs(values[1L]))
}

.empirical_variance_detect_constant_rates <- function(lambda_fun, mu_fun, grid, tol = 1e-10){
  lambda_vals <- .birth_death_eval_rate_fun(lambda_fun, grid)
  mu_vals <- .birth_death_eval_rate_fun(mu_fun, grid)

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
    return(.birth_death_li2_series_0_05(z, tol = tol, max_iter = max_iter))
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

    rec_t <- .birth_death_solve_reconstructed_process(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      t = t,
      n_steps = n_steps
    )

    phi_t <- tail(rec_t$phi, 1L)
    f_vals <- .empirical_variance_fprob_from_phi(rec_t$phi, phi_t, tol = tol)
    int_f <- .birth_death_trapz_vec(rec_t$grid, f_vals)

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


.empirical_variance_sim_BD_BM_var_path_inhom <- function(lambda_fun, mu_fun, sigma2,
                                     tmax, dt,
                                     envelope = NULL,
                                     n_envelope = 4000,
                                     safety_factor = 1.05){
  out <- .birth_death_brownian_sim_path_inhom(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    sigma2 = sigma2,
    tmax = tmax,
    dt = dt,
    x0 = 0,
    envelope = envelope,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  list(
    grid = out$grid,
    var = out$var,
    n = out$n
  )
}

.empirical_variance_simulate_BD_BM_paths <- function(lambda_fun, mu_fun, sigma2,
                                 tmax, dt, B,
                                 seed = NULL,
                                 n_envelope = 4000,
                                 safety_factor = 1.05){
  out <- .birth_death_brownian_simulate_paths(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    sigma2 = sigma2,
    tmax = tmax,
    dt = dt,
    B = B,
    x0 = 0,
    seed = seed,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  list(
    grid = out$grid,
    sigma2 = out$sigma2,
    V = out$V,
    N = out$N,
    B = out$B,
    tmax = out$tmax,
    dt = out$dt
  )
}

## 5) Empirical summaries extracted from simulations
## ================================================================

.empirical_variance_compute_empirical_BD_BM <- function(sim_res){
  summary <- .birth_death_brownian_summarise_replicates(sim_res)

  list(
    grid = summary$time,
    p_surv_emp = summary$p_survival_empirical,
    V_emp_uncond = summary$empirical_variance_expectation_empirical,
    V_emp_cond_same_time = summary$empirical_variance_expectation_empirical_cond_survival
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
