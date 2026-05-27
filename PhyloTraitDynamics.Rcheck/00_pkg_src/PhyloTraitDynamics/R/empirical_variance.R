## ================================================================
## Uses shared internal birth-death utilities.
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



#' Plot Expected Empirical Variance
#'
#' Plots simulated empirical variance paths, Monte Carlo summaries, and/or the
#' theoretical expectation on a common set of axes.
#'
#' @param simulation Optional `"empirical_variance_simulation"` object.
#' @param theory Optional data frame returned by
#'   `empirical_variance_compute_expectation()`.
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
                                                show_paths = TRUE,
                                                conditioning = c("none", "survival"),
                                                max_paths = Inf,
                                                alpha = 0.05,
                                                main = NULL,
                                                xlab = "time",
                                                ylab = "Empirical variance") {
  conditioning <- match.arg(conditioning)
  summary = NULL
  if (!is.null(simulation) && !inherits(simulation, "empirical_variance_simulation")) {
    stop("'simulation' must be NULL or the output of birth_death_brownian_simulate().", call. = FALSE)
  }

  if (is.null(simulation) && is.null(theory)) {
    stop("At least one of 'simulation' and 'theory' must be supplied.", call. = FALSE)
  }

  if (!is.null(simulation)) {
    summary <- simulation$summary
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


## ================================================================
## 0) Generic utilities
## ================================================================

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
      (2 * (1 + lambda * tt) / lambda) * .dilogarithm_li2_real(z, tol = li2_tol)

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
        (2 / lambda) * (1 - e_yule + e_yule * .dilogarithm_li2_real(z, tol = li2_tol))
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
    (.dilogarithm_li2_real(li2_arg_1, tol = li2_tol) -
       .dilogarithm_li2_real(li2_arg_2, tol = li2_tol))

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
## 1) Stable integrands for the theoretical expectations
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
## 2) Theory on a time grid
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
