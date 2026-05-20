## ================================================================
## Uses shared birth-death utilities from birth_death_utils.R.
## ================================================================

## ================================================================
## Public API for the variance of the empirical mean under a
## birth-death process with Brownian traits.
##
## This file is a package-oriented wrapper around the original script.
## The numerical and simulation routines are kept as internal functions
## with an .empirical_mean_ prefix.
## ================================================================

#' Variance of the Empirical Mean Under Birth-Death Brownian Dynamics
#'
#' Computes the theoretical variance of the empirical mean of Brownian traits,
#' conditional on survival at each time.
#'
#' @param birth Non-negative birth rate, either a numeric constant or a function
#'   of time.
#' @param death Non-negative death rate, either a numeric constant or a function
#'   of time.
#' @param sigma2 Brownian variance parameter.
#' @param time_end,time_step Time grid definition.
#' @param n_steps Number of integration steps for numerical theory.
#' @param tol_sing,li2_tol Numerical tolerances.
#'
#' @return A data frame with time, survival probability, empirical mean
#'   variance, and conditioning columns. The original script result is attached
#'   as attribute `"original_result"`.
#'
#' @export
empirical_mean_compute_variance <- function(birth,
                                            death,
                                            sigma2 = 1,
                                            time_end,
                                            time_step,
                                            n_steps = 1200,
                                            tol_sing = 1e-7,
                                            li2_tol = 1e-12) {
    conditioning = "survival"
	time_start = 0
  if (!identical(conditioning, "survival")) {
    stop("Only conditioning = 'survival' is implemented by the original script.", call. = FALSE)
  }

  grid <- .empirical_mean_make_public_grid(time_start, time_end, time_step)
  birth_fun <- .birth_death_make_rate_function(birth, "birth")
  death_fun <- .birth_death_make_rate_function(death, "death")

  res <- .empirical_mean_compute_theory_BD_BM_meanvar(
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


#' Summarise an Empirical Mean Simulation
#'
#' Computes Monte Carlo summaries from a simulation returned by
#' `birth_death_brownian_simulate()`.
#'
#' @param x An `"empirical_mean_simulation"` object.
#'
#' @return A data frame with time, empirical survival probability,
#'   survival-conditioned empirical mean, and survival-conditioned variance of
#'   the empirical mean.
#'
#' @keywords internal
empirical_mean_summarise_simulation <- function(x) {
  if (!inherits(x, "empirical_mean_simulation")) {
    stop("'x' must be the output of birth_death_brownian_simulate().", call. = FALSE)
  }

  if (!is.null(x$summary)) {
    return(x$summary)
  }

  .birth_death_brownian_summarise_replicates(x)
}


#' Plot Variance of the Empirical Mean
#'
#' Plots Monte Carlo summaries and/or theoretical values for the variance of the
#' empirical mean.
#'
#' @param simulation Optional `"empirical_mean_simulation"` object.
#' @param theory Optional data frame returned by
#'   `empirical_mean_compute_variance()`.
#' @param summary Optional data frame returned by
#'   `empirical_mean_summarise_simulation()`. If omitted and `simulation` is
#'   supplied, it is computed automatically.
#' @param main,xlab,ylab Base graphics labels.
#'
#' @return Invisibly returns a list containing `simulation`, `theory`, and
#'   `summary`.
#'
#' @export
empirical_mean_plot_variance <- function(simulation = NULL,
                                         theory = NULL,
                                         summary = NULL,
                                         main = NULL,
                                         xlab = "time",
                                         ylab = "Variance of empirical mean") {
  if (!is.null(simulation) && !inherits(simulation, "empirical_mean_simulation")) {
    stop("'simulation' must be NULL or the output of birth_death_brownian_simulate().", call. = FALSE)
  }

  if (is.null(summary) && !is.null(simulation)) {
    summary <- empirical_mean_summarise_simulation(simulation)
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

.empirical_mean_make_public_grid <- function(time_start = 0, time_end, time_step) {
  if (!is.numeric(time_start) || length(time_start) != 1L ||
      !is.finite(time_start) || time_start < 0) {
    stop("time_start must be a finite nonnegative number.")
  }
  if (!is.numeric(time_end) || length(time_end) != 1L ||
      !is.finite(time_end) || time_end < 0) {
    stop("time_end must be a finite nonnegative number.")
  }
  if (!is.numeric(time_step) || length(time_step) != 1L ||
      !is.finite(time_step) || time_step <= 0) {
    stop("time_step must be a strictly positive number.")
  }
  if (time_start > time_end) {
    stop("time_start must be smaller than or equal to time_end.")
  }

  grid <- seq(time_start, time_end, by = time_step)
  if (length(grid) == 0L || abs(tail(grid, 1L) - time_end) > 1e-14) {
    grid <- c(grid, time_end)
  }
  unique(grid)
}


## ================================================================
## Time-dependent birth-death process + Brownian traits
## Variance of the empirical mean, conditional on survival
## ---------------------------------------------------------------
## This script provides:
##   - numerical theory for Var( empirical mean at time t | N(t) > 0 ),
##   - inhomogeneous birth-death simulation by thinning,
##   - Brownian trait simulation on the living lineages,
##   - Monte Carlo estimation of Var( empirical mean | N(t) > 0 ),
##   - comparison and plotting utilities.
##
## There is deliberately:
##   - no closed-form special handling for constant rates,
##   - no conditioning on N(t) >= 2,
##   - no conditioning on survival at a final time T.
##
## A single simulated replicate produces a trajectory of the empirical mean,
## not a trajectory of its variance. The variance curve is estimated across
## Monte Carlo replicates, at each time point, conditional on N(t) > 0.
## ================================================================


## ================================================================
## 0) Generic utilities
## ================================================================

## ================================================================
## 1) Reconstructed process at fixed final time t
##    Solve backward:
##      g_t'(w) = (lambda(w) - mu(w)) g_t(w) - lambda(w),
##      g_t(t)  = 1.
## ================================================================

## ================================================================
## 2) Integrand for Var(empirical mean | N(t)>0)
##
## For a fixed target time t, let phi_s = phi_t(s) and phi_T = phi_t(t).
## The integrand is
##
## 1 + phi_T/(exp(phi_T)-1)
##   - 2 phi_T/(exp(phi_s)-1)
##   - 2 (exp(phi_T)-exp(phi_s))/(exp(phi_s)-1)^2
##       * [ Li_2(1-exp(-(phi_T-phi_s))) - Li_2(1-exp(-phi_T)) ].
##
## It has a removable singularity at phi_s = 0. Its limit there is 0.
## ================================================================

.empirical_mean_mean_variance_integrand_from_phi <- function(phi_s, phi_t,
                                             tol_sing = 1e-7,
                                             li2_tol = 1e-12){
  x <- as.numeric(phi_s)
  a <- as.numeric(phi_t)

  if (length(a) != 1L || !is.finite(a) || a < 0){
    stop("'phi_t' must be a single finite nonnegative number.")
  }
  if (any(!is.finite(x))){
    stop("Non-finite values in 'phi_s'.")
  }

  x <- pmin(pmax(x, 0), a)
  out <- numeric(length(x))

  ## If no reconstructed birth can occur, the integrand is zero and
  ## the variance is simply sigma2 * t.
  if (a < tol_sing){
    return(out)
  }

  keep <- x >= tol_sing
  if (!any(keep)){
    return(out)
  }

  xx <- x[keep]
  exp_a <- exp(a)
  exp_x <- exp(xx)
  denom_x <- expm1(xx)
  denom_a <- expm1(a)

  arg_delta <- 1 - exp(-pmax(a - xx, 0))
  arg_a <- 1 - exp(-a)

  dli <- .dilogarithm_li2_01(arg_delta, tol = li2_tol) -
    .dilogarithm_li2_01(arg_a, tol = li2_tol)

  out[keep] <-
    1 + a / denom_a -
    2 * a / denom_x -
    2 * (exp_a - exp_x) / (denom_x^2) * dli

  ## Very small negative values may occur by cancellation near the origin.
  out
}


## ================================================================
## 3) Theory on a time grid
## ================================================================

.empirical_mean_compute_theory_BD_BM_meanvar <- function(lambda_fun, mu_fun, sigma2, grid,
                                         n_steps = 1200,
                                         tol_sing = 1e-7,
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

  n_grid <- length(grid)
  p_surv <- numeric(n_grid)
  Var_mean_cond_surv <- numeric(n_grid)

  reconstructed <- if (return_reconstructed) vector("list", n_grid) else NULL

  for (i in seq_along(grid)){
    t <- grid[i]

    if (t == 0){
      p_surv[i] <- 1
      Var_mean_cond_surv[i] <- 0
      if (return_reconstructed){
        reconstructed[[i]] <- .birth_death_solve_reconstructed_process(
          lambda_fun = lambda_fun,
          mu_fun = mu_fun,
          t = 0,
          n_steps = 1
        )
      }
      next
    }

    rec_t <- .birth_death_solve_reconstructed_process(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      t = t,
      n_steps = n_steps
    )

    phi_t <- tail(rec_t$phi, 1L)
    integrand <- .empirical_mean_mean_variance_integrand_from_phi(
      phi_s = rec_t$phi,
      phi_t = phi_t,
      tol_sing = tol_sing,
      li2_tol = li2_tol
    )
    int_val <- .birth_death_trapz_vec(rec_t$grid, integrand)

    p_surv[i] <- rec_t$p_surv[1L]
    Var_mean_cond_surv[i] <- sigma2 * (t - int_val)

    ## Numerical cancellation may produce tiny negative values for very small t.
    if (Var_mean_cond_surv[i] < 0 && Var_mean_cond_surv[i] > -1e-10){
      Var_mean_cond_surv[i] <- 0
    }

    if (return_reconstructed){
      reconstructed[[i]] <- rec_t
    }
  }

  list(
    grid = grid,
    sigma2 = sigma2,
    p_surv = p_surv,
    Var_mean_th_cond_surv = Var_mean_cond_surv,
    reconstructed = reconstructed
  )
}


## ================================================================
## 4) Inhomogeneous birth-death Brownian simulation
##
## Simulation is now handled by birth_death_brownian_simulate() and
## its internal helpers in birth_death.R.
