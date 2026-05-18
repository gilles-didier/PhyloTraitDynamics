## ================================================================
## Shared internal utilities for birth-death computations.
##
## These functions are intentionally non-exported. They collect only
## low-level utilities shared by the MRCA-age, empirical-mean, and
## empirical-variance modules.
## ================================================================
.birth_death_make_rate_function <- function(rate, name = "rate") {
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

.birth_death_make_time_grid <- function(tmax, dt){
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

.birth_death_trapz_vec <- function(x, y){
  n <- length(x)
  if (n != length(y)){
    stop("'x' and 'y' must have the same length.")
  }
  if (n <= 1L){
    return(0)
  }
  sum(diff(x) * (y[-n] + y[-1L]) / 2)
}

.birth_death_cumtrapz_vec <- function(x, y){
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

.birth_death_eval_rate_fun <- function(rate_fun, x){
  out <- rate_fun(x)
  if (length(out) == 1L && length(x) > 1L){
    out <- rep(out, length(x))
  }
  if (length(out) != length(x)){
    stop("A rate function must return either a scalar or a vector of the same length as its input.")
  }
  as.numeric(out)
}

.birth_death_check_rate_functions <- function(lambda_fun, mu_fun, grid){
  lambda_vals <- .birth_death_eval_rate_fun(lambda_fun, grid)
  mu_vals <- .birth_death_eval_rate_fun(mu_fun, grid)

  if (any(!is.finite(lambda_vals)) || any(!is.finite(mu_vals))){
    stop("Non-finite rate values detected.")
  }
  if (any(lambda_vals < 0) || any(mu_vals < 0)){
    stop("Rates must be nonnegative on the working grid.")
  }

  invisible(list(lambda = lambda_vals, mu = mu_vals))
}

.birth_death_poly_rate_fun <- function(coeffs, floor_value = 0){
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

.birth_death_constant_rate_fun <- function(rate){
  if (!is.finite(rate) || rate < 0){
    stop("'rate' must be a finite nonnegative number.")
  }
  force(rate)
  function(t){
    rep(rate, length(as.numeric(t)))
  }
}

.birth_death_li2_series_0_05 <- function(x, tol = 1e-12, max_iter = 10000L){
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

.birth_death_solve_reconstructed_process <- function(lambda_fun, mu_fun, t,
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
    lambda0 <- .birth_death_eval_rate_fun(lambda_fun, 0)
    mu0 <- .birth_death_eval_rate_fun(mu_fun, 0)

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
    .birth_death_check_rate_functions(lambda_fun, mu_fun, grid_desc)
  }

  rhs_g <- function(w, g){
    (.birth_death_eval_rate_fun(lambda_fun, w) - .birth_death_eval_rate_fun(mu_fun, w)) * g -
      .birth_death_eval_rate_fun(lambda_fun, w)
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
  lambda_vals <- .birth_death_eval_rate_fun(lambda_fun, grid)
  mu_vals <- .birth_death_eval_rate_fun(mu_fun, grid)

  rebirth <- lambda_vals / g
  phi <- .birth_death_cumtrapz_vec(grid, rebirth)

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

.birth_death_make_rate_envelope <- function(lambda_fun, mu_fun, tmax,
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

  checked <- .birth_death_check_rate_functions(lambda_fun, mu_fun, grid)
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

.birth_death_sample_next_event_thinning <- function(time, k, tmax,
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

      lambda_now <- .birth_death_eval_rate_fun(lambda_fun, current)
      mu_now <- .birth_death_eval_rate_fun(mu_fun, current)
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

## ================================================================
## Shared birth-death Brownian simulation utilities.
##
## These functions simulate the birth-death tree and Brownian traits
## once, and record all empirical quantities needed by the empirical-
## mean and empirical-variance modules.
## ================================================================
.birth_death_brownian_sim_path_inhom <- function(lambda_fun, mu_fun, sigma2,
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

  grid <- .birth_death_make_time_grid(tmax, dt)
  n_grid <- length(grid)

  if (is.null(envelope)){
    envelope <- .birth_death_make_rate_envelope(
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
  var_path <- numeric(n_grid)
  n_path <- integer(n_grid)

  mean_path[1L] <- x0
  var_path[1L] <- 0
  n_path[1L] <- 1L

  next_event <- .birth_death_sample_next_event_thinning(
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
      var_path[idx] <- if (k >= 2L) stats::var(traits) else 0
      next
    }

    if (abs(time - next_event) < 1e-12){
      lambda_now <- .birth_death_eval_rate_fun(lambda_fun, time)
      mu_now <- .birth_death_eval_rate_fun(mu_fun, time)
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

      next_event <- .birth_death_sample_next_event_thinning(
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
    var = var_path,
    n = n_path
  )
}

.birth_death_brownian_simulate_paths <- function(lambda_fun, mu_fun, sigma2,
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

  grid <- .birth_death_make_time_grid(tmax, dt)
  n_grid <- length(grid)

  envelope <- .birth_death_make_rate_envelope(
    lambda_fun = lambda_fun,
    mu_fun = mu_fun,
    tmax = tmax,
    n_envelope = n_envelope,
    safety_factor = safety_factor
  )

  M <- matrix(NA_real_, nrow = B, ncol = n_grid)
  V <- matrix(0, nrow = B, ncol = n_grid)
  N <- matrix(0L, nrow = B, ncol = n_grid)

  for (b in seq_len(B)){
    out <- .birth_death_brownian_sim_path_inhom(
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
    V[b, ] <- out$var
    N[b, ] <- out$n
  }

  out <- list(
    grid = grid,
    sigma2 = sigma2,
    x0 = x0,
    M = M,
    V = V,
    N = N,
    B = B,
    tmax = tmax,
    dt = dt
  )
  out$summary <- .birth_death_brownian_summarise_replicates(out)
  out
}

.birth_death_brownian_summarise_replicates <- function(sim_res){
  grid <- sim_res$grid
  N <- sim_res$N

  if (is.null(grid) || is.null(N)){
    stop("'sim_res' must contain at least 'grid' and 'N'.")
  }

  has_M <- !is.null(sim_res$M)
  has_V <- !is.null(sim_res$V)

  M <- if (has_M) sim_res$M else matrix(NA_real_, nrow = nrow(N), ncol = ncol(N))
  V <- if (has_V) sim_res$V else matrix(NA_real_, nrow = nrow(N), ncol = ncol(N))

  alive_same_time <- N > 0
  p_surv_emp <- colMeans(alive_same_time)

  Mean_emp_cond_surv <- rep(NA_real_, length(grid))
  Var_mean_emp_cond_surv <- rep(NA_real_, length(grid))
  V_emp_cond_same_time <- rep(NA_real_, length(grid))

  V_emp_uncond <- if (has_V) colMeans(V) else rep(NA_real_, length(grid))

  for (j in seq_along(grid)){
    keep <- alive_same_time[, j]
    n_keep <- sum(keep)

    if (has_M && n_keep > 0L){
      m <- M[keep, j]
      Mean_emp_cond_surv[j] <- mean(m)
    }

    if (has_M && n_keep >= 2L){
      Var_mean_emp_cond_surv[j] <- stats::var(M[keep, j])
    }

    if (has_V && n_keep > 0L){
      V_emp_cond_same_time[j] <- mean(V[keep, j])
    }
  }

  data.frame(
    time = grid,
    p_survival_empirical = p_surv_emp,
    empirical_mean_empirical_cond_survival = Mean_emp_cond_surv,
    empirical_mean_variance_empirical_cond_survival = Var_mean_emp_cond_surv,
    empirical_variance_expectation_empirical = V_emp_uncond,
    empirical_variance_expectation_empirical_cond_survival = V_emp_cond_same_time
  )
}

