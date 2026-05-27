## ============================================================
## Uses shared internal birth-death utilities.
## ============================================================

## ============================================================
## MRCA age distribution through time under generalized BD rates
## ============================================================

#' MRCA Age Distribution at One Time
#'
#' Computes the distribution and summary statistics of the most recent common
#' ancestor age at a fixed time under a birth-death process.
#'
#' @param birth Non-negative birth rate, either a numeric constant or a function
#'   of time.
#' @param death Non-negative death rate, either a numeric constant or a function
#'   of time.
#' @param time Observation time.
#' @param time_step Internal time step. If `NULL`, a default step is derived
#'   from `time`.
#' @param probs Quantile probabilities to report.
#' @param n_internal Optional number of internal subintervals; when supplied it
#'   overrides `time_step`.
#'
#' @return An object of class `"mrca_age_distribution"` with distribution grid,
#'   summary statistics, requested probabilities, and the original result.
#'
#' @export
mrca_age_compute_distribution <- function(birth,
                                          death,
                                          time,
                                          time_step = NULL,
                                          probs = c(0.25, 0.5, 0.75),
                                          n_internal = NULL) {
  if (!is.numeric(time) || length(time) != 1L) {
    stop("time must be a number.")
  }

  if (!is.null(n_internal)) {
    if (!is.numeric(n_internal) || length(n_internal) != 1L || n_internal <= 0) {
      stop("n_internal must be a positive number.")
    }
    time_step <- time / n_internal
  }

  if (time > 0) {
    if (is.null(time_step)) {
      time_step <- time / 50
    }

    if (!is.numeric(time_step) || length(time_step) != 1L || time_step <= 0) {
      stop("time_step must be a positive number.")
    }
  } else {
    if (is.null(time_step)) {
      time_step <- 1
    }
  }

  if (any(probs <= 0 | probs >= 1)) {
    stop("probs must lie strictly between 0 and 1.")
  }

  probs <- sort(probs)

  birth_fun <- .birth_death_make_rate_function(birth, "birth")
  death_fun <- .birth_death_make_rate_function(death, "death")

  out <- .mrca_age_distribution_one_time(
    time,
    birth_fun = birth_fun,
    death_fun = death_fun,
    dt_internal = time_step,
    probs = probs
  )

  result <- list(
    summary = as.data.frame(as.list(out$summary)),
    grid = out$grid,
    probs = probs,
    original_result = out,
    call = match.call()
  )

  class(result) <- "mrca_age_distribution"

  result
}


#' MRCA Age Distribution Dynamics
#'
#' Computes MRCA age distribution summaries over a time grid under a birth-death
#' process.
#'
#' @param birth Non-negative birth rate, either a numeric constant or a function
#'   of time.
#' @param death Non-negative death rate, either a numeric constant or a function
#'   of time.
#' @param time_end,time_step Time grid definition.
#' @param probs Quantile probabilities to report.
#' @param n_internal_per_step Number of internal integration substeps per output
#'   time step.
#' @param include_zero Logical; whether to include time zero in the original
#'   dynamics computation.
#'
#' @return An object of class `"mrca_age_dynamic"` containing the summary data
#'   frame, distribution grids, probabilities, and original call.
#'
#' @export
mrca_age_compute_dynamics <- function(birth,
                                      death,
                                      time_end,
                                      time_step,
                                      probs = c(0.25, 0.50, 0.75),
                                      n_internal_per_step = 10L,
                                      include_zero = TRUE) {
  if (!is.numeric(time_end) || length(time_end) != 1L || time_end <= 0) {
    stop("time_end must be a positive number.")
  }

  if (!is.numeric(time_step) || length(time_step) != 1L || time_step <= 0) {
    stop("time_step must be a positive number.")
  }

  if (any(probs <= 0 | probs >= 1)) {
    stop("probs must lie strictly between 0 and 1.")
  }

  probs <- sort(probs)

  birth_fun <- .birth_death_make_rate_function(birth, "birth")
  death_fun <- .birth_death_make_rate_function(death, "death")

  dt_internal <- time_step / n_internal_per_step

  times <- seq(time_step, time_end, by = time_step)

  if (include_zero) {
    times <- c(0, times)
  }

  if (tail(times, 1L) < time_end) {
    times <- c(times, time_end)
  }

  out <- lapply(
    times,
    .mrca_age_distribution_one_time,
    birth_fun = birth_fun,
    death_fun = death_fun,
    dt_internal = dt_internal,
    probs = probs
  )

  summary <- as.data.frame(
    do.call(rbind, lapply(out, `[[`, "summary"))
  )

  rownames(summary) <- NULL

  grids <- lapply(seq_along(out), function(i) {
    cbind(time = times[i], out[[i]]$grid)
  })

  result <- list(
    summary = summary,
    grids = grids,
    probs = probs,
    call = match.call()
  )

  class(result) <- "mrca_age_dynamic"

  result
}


#' Plot MRCA Age Distribution Dynamics
#'
#' Plots MRCA age summaries over time, including median, optional interquartile
#' band, and optional expectation.
#'
#' @param x An `"mrca_age_dynamic"` object.
#' @param show_iqr Logical; whether to draw the interquartile band.
#' @param show_expectation Logical; whether to draw the expectation curve.
#' @param main,xlab,ylab Base graphics labels.
#' @param add_legend Logical; whether to draw the plot legend.
#' @param ... Additional graphical arguments.
#'
#' @return Invisibly returns `x`.
#'
#' @export
mrca_age_plot_dynamics <- function(x,
                                   show_iqr = TRUE,
                                   show_expectation = TRUE,
                                   main = NA,
                                   xlab = "Time",
                                   ylab = "MRCA age",
                                   add_legend = TRUE,
                                   ...) {
  if (!inherits(x, "mrca_age_dynamic")) {
    stop("x must be the output of mrca_age_compute_dynamics().")
  }

  d <- x$summary

  qnames <- paste0("q", x$probs * 100)

  if (!all(c("time", "mean") %in% names(d)) || !all(qnames %in% names(d))) {
    stop("The summary table has unexpected column names.")
  }

  ylim <- range(d[, c("mean", qnames)], na.rm = TRUE)

  if (!all(is.finite(ylim))) {
    ylim <- c(0, max(d$time, na.rm = TRUE))
  }

  plot(
    d$time,
    d$mean,
    type = "n",
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  q25 <- qnames[1L]
  q50 <- qnames[which.min(abs(x$probs - 0.5))]
  q75 <- qnames[length(qnames)]

  if (show_iqr && length(qnames) >= 3L) {
    polygon(
      c(d$time, rev(d$time)),
      c(d[[q25]], rev(d[[q75]])),
      col = rgb(1,0,0,alpha = 0.2),
      border = NA
    )
  }

  lines(d$time, d[[q50]], col = "red", lty = 1, lwd = 1.5)
  if (length(qnames) >= 3L) {
    lines(d$time, d[[q25]], col = "red", lty = 2, lwd = 1)
    lines(d$time, d[[q75]], col = "red", lty = 2, lwd = 1)
  }
  if (show_expectation) {
    lines(d$time, d$mean, col = "blue", lty = 1, lwd = 0.5)
  }
  if(add_legend) {
	  legend(
		"topleft",
		legend = c("Expectation", "Median", "Quartiles"),
		lty = c(1,1,2),
		lwd = c(0.5, 1.5, 1),
	#    pch = c(NA, NA, 15),
	#    pt.cex = c(NA, NA, 2),
		col = c("blue", "red", "red"),
		bty = "n"
	  )
  }

  invisible(x)
}


.mrca_age_cdf_from_theta <- function(theta, theta_T, eps = 1e-7) {
  if (!is.finite(theta_T) || theta_T <= 0) {
    return(rep(NA_real_, length(theta)))
  }

  x <- pmin(pmax(theta, 0), theta_T)

  ## Stable version of
  ##
  ## g(x) = (1 - exp(x) + x exp(x)) / (exp(x) - 1)^2
  ##
  ## with expansion near 0:
  ##
  ## g(x) = 1/2 - x/6 + x^3/180 + O(x^5)
  ##
  ## This avoids numerical instability at the origin.
  g <- numeric(length(x))
  small <- abs(x) < eps

  g[small] <- 0.5 - x[small] / 6 + x[small]^3 / 180

  ex <- exp(x[!small])
  g[!small] <- (1 - ex + x[!small] * ex) / expm1(x[!small])^2

  ratio <- (expm1(theta_T) - expm1(x)) / expm1(theta_T)

  F <- 1 - 2 * ratio * g

  F[x <= 0] <- 0
  F[x >= theta_T] <- 1

  F <- pmin(pmax(F, 0), 1)

  ## Small numerical errors can create tiny non-monotonicities.
  F <- cummax(F)
  F[length(F)] <- 1

  F
}


.mrca_age_quantiles_from_cdf_grid <- function(s, F, probs) {
  ok <- is.finite(s) & is.finite(F)

  s <- s[ok]
  F <- F[ok]

  if (length(s) < 2L || max(F) < min(probs)) {
    return(rep(NA_real_, length(probs)))
  }

  F <- cummax(pmin(pmax(F, 0), 1))

  keep <- !duplicated(F)

  approx(
    x = F[keep],
    y = s[keep],
    xout = probs,
    rule = 2
  )$y
}


.mrca_age_distribution_one_time <- function(T, birth_fun, death_fun, dt_internal, probs) {
  if (T <= 0) {
    qs <- setNames(rep(NA_real_, length(probs)), paste0("q", probs * 100))

    return(list(
      summary = c(
        time = T,
        mean = NA_real_,
        qs,
        theta_T = 0,
        p_survival_0_T = NA_real_
      ),
      grid = data.frame(
        s = 0,
        p_survival = NA_real_,
        beta = NA_real_,
        theta = 0,
        cdf = NA_real_
      )
    ))
  }

  n_grid <- max(50L, ceiling(T / dt_internal) + 1L)

  s <- seq(0, T, length.out = n_grid)

  birth <- birth_fun(s)
  death <- death_fun(s)

  if (any(!is.finite(birth)) || any(!is.finite(death))) {
    stop("Rates must be finite on the time grid.")
  }

  if (any(birth < 0) || any(death < 0)) {
    stop("Birth and death rates must be non-negative.")
  }

  ## A(s) = integral_0^s (birth(u) - death(u)) du
  A <- .birth_death_cumtrapz_vec(s, birth - death)

  ## Survival probability:
  ##
  ## p_survival(s,T)
  ## = 1 / [1 + integral_s^T exp(-(A(u)-A(s))) death(u) du]
  ##
  ## Write the integral as:
  ##
  ## exp(A(s)) * integral_s^T exp(-A(u)) death(u) du
  h <- exp(-A) * death
  H <- .birth_death_cumtrapz_vec(s, h)

  tail_integral <- tail(H, 1L) - H

  p_survival <- 1 / (1 + exp(A) * tail_integral)

  ## Reconstructed birth rate beta_T(s)
  beta <- birth * p_survival

  ## theta(s) = integral_0^s beta_T(w) dw
  theta <- .birth_death_cumtrapz_vec(s, beta)
  theta_T <- tail(theta, 1L)

  if (!is.finite(theta_T) || theta_T <= 0) {
    F <- rep(NA_real_, length(s))
    mean_mrca <- NA_real_
    qs <- rep(NA_real_, length(probs))
  } else {
    F <- .mrca_age_cdf_from_theta(theta, theta_T)

    ## For a nonnegative variable supported on [0,T]:
    ## E[U_T] = integral_0^T (1 - F_T(s)) ds
    mean_mrca <- .birth_death_trapz_vec(s, 1 - F)

    qs <- .mrca_age_quantiles_from_cdf_grid(s, F, probs)
  }

  names(qs) <- paste0("q", probs * 100)

  list(
    summary = c(
      time = T,
      mean = mean_mrca,
      qs,
      theta_T = theta_T,
      p_survival_0_T = p_survival[1L]
    ),
    grid = data.frame(
      s = s,
      p_survival = p_survival,
      beta = beta,
      theta = theta,
      cdf = F
    )
  )
}
