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

make_time_grid <- function(tmax, dt){
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

trapz_vec <- function(x, y){
  n <- length(x)
  if (n != length(y)){
    stop("'x' and 'y' must have the same length.")
  }
  if (n <= 1L){
    return(0)
  }
  sum(diff(x) * (y[-n] + y[-1L]) / 2)
}

cumtrapz_vec <- function(x, y){
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

eval_rate_fun <- function(rate_fun, x){
  out <- rate_fun(x)
  if (length(out) == 1L && length(x) > 1L){
    out <- rep(out, length(x))
  }
  if (length(out) != length(x)){
    stop("A rate function must return either a scalar or a vector of the same length as its input.")
  }
  as.numeric(out)
}

check_rate_functions <- function(lambda_fun, mu_fun, grid){
  lambda_vals <- eval_rate_fun(lambda_fun, grid)
  mu_vals <- eval_rate_fun(mu_fun, grid)

  if (any(!is.finite(lambda_vals)) || any(!is.finite(mu_vals))){
    stop("Non-finite rate values detected.")
  }
  if (any(lambda_vals < 0) || any(mu_vals < 0)){
    stop("Rates must be nonnegative on the working grid.")
  }

  invisible(list(lambda = lambda_vals, mu = mu_vals))
}

poly_rate_fun <- function(coeffs, floor_value = 0){
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

constant_rate_fun <- function(rate){
  if (!is.finite(rate) || rate < 0){
    stop("'rate' must be a finite nonnegative number.")
  }
  force(rate)
  function(t){
    rep(rate, length(as.numeric(t)))
  }
}


## ================================================================
## 0b) Real dilogarithm Li_2(z) on 0 <= z <= 1
##     This is the only domain needed by the integral below.
## ================================================================

li2_series_0_05 <- function(x, tol = 1e-12, max_iter = 10000L){
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

li2_01_scalar <- function(z, tol = 1e-12, max_iter = 10000L){
  if (!is.finite(z)){
    stop("Non-finite argument passed to li2_01().")
  }
  if (z < -1e-12 || z > 1 + 1e-12){
    stop("li2_01() is implemented only for 0 <= z <= 1.")
  }

  z <- min(1, max(0, z))

  if (z == 0){
    return(0)
  }
  if (z == 1){
    return(pi^2 / 6)
  }
  if (z <= 0.5){
    return(li2_series_0_05(z, tol = tol, max_iter = max_iter))
  }

  ## Li_2(z) = pi^2/6 - log(z) log(1-z) - Li_2(1-z)
  pi^2 / 6 - log(z) * log1p(-z) -
    li2_series_0_05(1 - z, tol = tol, max_iter = max_iter)
}

li2_01 <- function(z, tol = 1e-12, max_iter = 10000L){
  vapply(
    as.numeric(z),
    li2_01_scalar,
    numeric(1),
    tol = tol,
    max_iter = as.integer(max_iter)
  )
}


## ================================================================
## 1) Reconstructed process at fixed final time t
##    Solve backward:
##      g_t'(w) = (lambda(w) - mu(w)) g_t(w) - lambda(w),
##      g_t(t)  = 1.
## ================================================================

solve_reconstructed_process <- function(lambda_fun, mu_fun, t,
                                        n_steps = 1200,
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
    lambda0 <- eval_rate_fun(lambda_fun, 0)
    mu0 <- eval_rate_fun(mu_fun, 0)

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
    check_rate_functions(lambda_fun, mu_fun, grid_desc)
  }

  rhs_g <- function(w, g){
    (eval_rate_fun(lambda_fun, w) - eval_rate_fun(mu_fun, w)) * g -
      eval_rate_fun(lambda_fun, w)
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
  lambda_vals <- eval_rate_fun(lambda_fun, grid)
  mu_vals <- eval_rate_fun(mu_fun, grid)

  rebirth <- lambda_vals / g
  phi <- cumtrapz_vec(grid, rebirth)

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

mean_variance_integrand_from_phi <- function(phi_s, phi_t,
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

  dli <- li2_01(arg_delta, tol = li2_tol) -
    li2_01(arg_a, tol = li2_tol)

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

compute_theory_BD_BM_meanvar <- function(lambda_fun, mu_fun, sigma2, grid,
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
        reconstructed[[i]] <- solve_reconstructed_process(
          lambda_fun = lambda_fun,
          mu_fun = mu_fun,
          t = 0,
          n_steps = 1
        )
      }
      next
    }

    rec_t <- solve_reconstructed_process(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      t = t,
      n_steps = n_steps
    )

    phi_t <- tail(rec_t$phi, 1L)
    integrand <- mean_variance_integrand_from_phi(
      phi_s = rec_t$phi,
      phi_t = phi_t,
      tol_sing = tol_sing,
      li2_tol = li2_tol
    )
    int_val <- trapz_vec(rec_t$grid, integrand)

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
## 4) Inhomogeneous birth-death simulation by thinning
## ================================================================

make_rate_envelope <- function(lambda_fun, mu_fun, tmax,
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

  checked <- check_rate_functions(lambda_fun, mu_fun, grid)
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

sample_next_event_thinning <- function(time, k, tmax,
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

      lambda_now <- eval_rate_fun(lambda_fun, current)
      mu_now <- eval_rate_fun(mu_fun, current)
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

sim_BD_BM_mean_path_inhom <- function(lambda_fun, mu_fun, sigma2,
                                      tmax, dt,
                                      x0 = 0,
                                      envelope = NULL,
                                      n_envelope = 4000,
                                      safety_factor = 1.05){
  if (!is.finite(sigma2) || sigma2 < 0){
    stop("'sigma2' must be a finite nonnegative number.")
  }
  if (!is.finite(x0)){
    stop("'x0' must be finite.")
  }

  grid <- make_time_grid(tmax, dt)
  n_grid <- length(grid)

  if (is.null(envelope)){
    envelope <- make_rate_envelope(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      tmax = tmax,
      n_envelope = n_envelope,
      safety_factor = safety_factor
    )
  }

  time <- 0
  traits <- x0
  k <- 1L
  idx <- 1L

  mean_path <- rep(NA_real_, n_grid)
  n_path <- integer(n_grid)
  mean_path[1L] <- x0
  n_path[1L] <- 1L

  next_event <- sample_next_event_thinning(
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
      mean_path[idx] <- if (k >= 1L) mean(traits) else NA_real_
      next
    }

    if (abs(time - next_event) < 1e-12){
      lambda_now <- eval_rate_fun(lambda_fun, time)
      mu_now <- eval_rate_fun(mu_fun, time)
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

      next_event <- sample_next_event_thinning(
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
    mean = mean_path,
    n = n_path
  )
}

simulate_BD_BM_mean_paths <- function(lambda_fun, mu_fun, sigma2,
                                      tmax, dt, B,
                                      x0 = 0,
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

  grid <- make_time_grid(tmax, dt)
  n_grid <- length(grid)

  envelope <- make_rate_envelope(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    tmax = tmax,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  M <- matrix(NA_real_, nrow = B, ncol = n_grid)
  N <- matrix(0L, nrow = B, ncol = n_grid)

  for (b in seq_len(B)){
    out <- sim_BD_BM_mean_path_inhom(
      lambda_fun = lambda_fun,
      mu_fun = mu_fun,
      sigma2 = sigma2,
      tmax = tmax,
      dt = dt,
      x0 = x0,
      envelope = envelope,
      safety_factor = safety_factor
    )
    M[b, ] <- out$mean
    N[b, ] <- out$n
  }

  list(
    grid = grid,
    sigma2 = sigma2,
    x0 = x0,
    M = M,
    N = N,
    B = B,
    tmax = tmax,
    dt = dt
  )
}


## ================================================================
## 5) Empirical summaries extracted from simulations
## ================================================================

compute_empirical_BD_BM_meanvar <- function(sim_res){
  grid <- sim_res$grid
  M <- sim_res$M
  N <- sim_res$N
  x0 <- if (!is.null(sim_res$x0)) sim_res$x0 else 0

  alive <- N > 0
  p_surv_emp <- colMeans(alive)

  ## Since E[bar X(t) | N(t)>0] = x0, the Monte Carlo estimator below uses
  ## the known centering x0 rather than the sample mean. This estimates
  ## E[(bar X(t)-x0)^2 | N(t)>0].
  Var_mean_emp_cond_surv <- rep(NA_real_, length(grid))
  Mean_emp_cond_surv <- rep(NA_real_, length(grid))

  for (j in seq_along(grid)){
    keep <- alive[, j]
    if (any(keep)){
      x <- M[keep, j]
      Var_mean_emp_cond_surv[j] <- mean((x - x0)^2)
      Mean_emp_cond_surv[j] <- mean(x)
    }
  }

  list(
    grid = grid,
    p_surv_emp = p_surv_emp,
    Mean_emp_cond_surv = Mean_emp_cond_surv,
    Var_mean_emp_cond_surv = Var_mean_emp_cond_surv
  )
}


## ================================================================
## 6) Comparison helpers
## ================================================================

compare_BD_BM_meanvar_fit <- function(sim_res, theory_res = NULL,
                                      empirical_res = NULL){
  if (is.null(empirical_res)){
    empirical_res <- compute_empirical_BD_BM_meanvar(sim_res = sim_res)
  }

  out <- data.frame(
    curve = c("variance_mean_cond_survival", "survival_probability"),
    max_abs_diff = NA_real_,
    rmse = NA_real_,
    final_emp = c(
      tail(empirical_res$Var_mean_emp_cond_surv, 1L),
      tail(empirical_res$p_surv_emp, 1L)
    ),
    final_th = NA_real_
  )

  if (is.null(theory_res)){
    return(out)
  }

  diff_var <- empirical_res$Var_mean_emp_cond_surv -
    theory_res$Var_mean_th_cond_surv
  diff_surv <- empirical_res$p_surv_emp - theory_res$p_surv

  out$max_abs_diff[1L] <- max(abs(diff_var), na.rm = TRUE)
  out$rmse[1L] <- sqrt(mean(diff_var^2, na.rm = TRUE))
  out$final_th[1L] <- tail(theory_res$Var_mean_th_cond_surv, 1L)

  out$max_abs_diff[2L] <- max(abs(diff_surv), na.rm = TRUE)
  out$rmse[2L] <- sqrt(mean(diff_surv^2, na.rm = TRUE))
  out$final_th[2L] <- tail(theory_res$p_surv, 1L)

  out
}


## ================================================================
## 7) Plotting
## ================================================================

plot_BD_BM_meanvar_results <- function(sim_res,
                                       theory_res = NULL,
                                       empirical_res = NULL,
                                       main = NULL,
                                       xlab = "Time",
                                       ylab = "Empirical mean"){
  grid <- sim_res$grid

  if (is.null(empirical_res)){
    empirical_res <- compute_empirical_BD_BM_meanvar(sim_res = sim_res)
  }

  y_candidates <- c(0, empirical_res$Var_mean_emp_cond_surv)
  if (!is.null(theory_res)){
    y_candidates <- c(y_candidates, theory_res$Var_mean_th_cond_surv)
  }

  ymax <- max(y_candidates, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0){
    ymax <- 1
  }

  plot(
    grid,
    empirical_res$Var_mean_emp_cond_surv,
    type = "l",
    lwd = 4,
    xlab = xlab,
    ylab = ylab,
    ylim = c(0, ymax),
    main = main
  )

  legend_labels <- "simulations variance"
  legend_lty <- 1
  legend_lwd <- 3
  legend_col <- "black"

  if (!is.null(theory_res)){
    lines(grid, theory_res$Var_mean_th_cond_surv,
          col = "green", lwd = 3, lty = 2)
    legend_labels <- c(legend_labels, "theoretical variance")
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 3)
    legend_col <- c(legend_col, "green")
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

plot_BD_BM_meanvar_survival <- function(sim_res,
                                        theory_res = NULL,
                                        empirical_res = NULL,
                                        main = NULL,
                                        xlab = "time",
                                        ylab = "Survival probability"){
  grid <- sim_res$grid

  if (is.null(empirical_res)){
    empirical_res <- compute_empirical_BD_BM_meanvar(sim_res = sim_res)
  }

  plot(
    grid,
    empirical_res$p_surv_emp,
    type = "l",
    lwd = 3,
    xlab = xlab,
    ylab = ylab,
    ylim = c(0, 1),
    main = main
  )

  legend_labels <- "Monte Carlo"
  legend_lty <- 1
  legend_lwd <- 3
  legend_col <- "black"

  if (!is.null(theory_res)){
    lines(grid, theory_res$p_surv, col = "forestgreen", lwd = 2, lty = 2)
    legend_labels <- c(legend_labels, "theory")
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 2)
    legend_col <- c(legend_col, "forestgreen")
  }

  legend(
    "topright",
    legend = legend_labels,
    lty = legend_lty,
    lwd = legend_lwd,
    col = legend_col,
    bty = "n"
  )

  invisible(list(sim = sim_res, empirical = empirical_res, theory = theory_res))
}

plot_BD_BM_empirical_mean_paths <- function(sim_res,
                                            max_paths = 100,
                                            alpha = 0.08,
                                            main = NULL,
                                            xlab = "time",
                                            ylab = "Empirical mean"){
  grid <- sim_res$grid
  M <- sim_res$M
  B <- nrow(M)
  path_idx <- seq_len(min(as.integer(max_paths), B))

  y_candidates <- as.numeric(M[path_idx, , drop = FALSE])
  yrange <- range(y_candidates, na.rm = TRUE)
  if (!all(is.finite(yrange)) || diff(yrange) == 0){
    yrange <- yrange + c(-1, 1)
  }

  plot(
    grid, rep(NA_real_, length(grid)),
    type = "n",
    xlab = xlab,
    ylab = ylab,
    ylim = yrange,
    main = main
  )

  col_path <- adjustcolor("black", alpha.f = alpha)
  for (b in path_idx){
    lines(grid, M[b, ], col = col_path)
  }

  invisible(sim_res)
}


## ================================================================
## 8) Convenience wrapper
## ================================================================

run_BD_BM_meanvar_experiment <- function(lambda_fun, mu_fun, sigma2,
                                         tmax, dt, B,
                                         x0 = 0,
                                         seed = NULL,
                                         n_steps_theory = 1200,
                                         n_envelope = 4000,
                                         safety_factor = 1.05,
                                         tol_sing = 1e-7,
                                         li2_tol = 1e-12){
  sim_res <- simulate_BD_BM_mean_paths(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    sigma2 = sigma2,
    tmax = tmax,
    dt = dt,
    B = B,
    x0 = x0,
    seed = seed,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  theory_res <- compute_theory_BD_BM_meanvar(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    sigma2 = sigma2,
    grid = sim_res$grid,
    n_steps = n_steps_theory,
    tol_sing = tol_sing,
    li2_tol = li2_tol
  )

  empirical_res <- compute_empirical_BD_BM_meanvar(sim_res = sim_res)

  list(
    sim = sim_res,
    theory = theory_res,
    empirical = empirical_res
  )
}


## ================================================================
## 9) Examples
## ================================================================

## Time-dependent example.
# lambda_fun <- poly_rate_fun(c(0.6, 0.8, -0.25))
# mu_fun <- poly_rate_fun(c(0.1, 0.05))
#
# res <- run_BD_BM_meanvar_experiment(
#   lambda_fun = lambda_fun,
#   mu_fun = mu_fun,
#   sigma2 = 1,
#   tmax = 3,
#   dt = 0.02,
#   B = 2000,
#   seed = 3,
#   n_steps_theory = 1200,
#   n_envelope = 4000
# )
#
# compare_BD_BM_meanvar_fit(res$sim, res$theory, res$empirical)
#
# plot_BD_BM_meanvar_results(
#   sim_res = res$sim,
#   theory_res = res$theory,
#   empirical_res = res$empirical,
#   main = "Variance of empirical mean | N(t)>0"
# )
#
# plot_BD_BM_meanvar_survival(
#   sim_res = res$sim,
#   theory_res = res$theory,
#   empirical_res = res$empirical,
#   main = "Survival probability"
# )

## Constant-rate example, still using the same numerical theory.
# lambda_fun <- constant_rate_fun(0.8)
# mu_fun <- constant_rate_fun(0.2)
#
# res_const <- run_BD_BM_meanvar_experiment(
#   lambda_fun = lambda_fun,
#   mu_fun = mu_fun,
#   sigma2 = 1,
#   tmax = 3,
#   dt = 0.02,
#   B = 2000,
#   seed = 4,
#   n_steps_theory = 1200,
#   n_envelope = 4000
# )
#
# compare_BD_BM_meanvar_fit(res_const$sim, res_const$theory, res_const$empirical)
#
# plot_BD_BM_meanvar_results(
#   sim_res = res_const$sim,
#   theory_res = res_const$theory,
#   empirical_res = res_const$empirical,
#   main = "Constant-rate BD: Var(empirical mean | N(t)>0)"
# )
