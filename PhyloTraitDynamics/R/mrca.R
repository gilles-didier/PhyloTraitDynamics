## Public API for MRCA age distributions.

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
#' @param time_step Internal time step. If `NULL`, a default step is derived from
#'   `time`.
#' @param probs Quantile probabilities to report.
#' @param n_internal Optional number of internal subintervals; when supplied it
#'   overrides `time_step`.
#'
#' @return An object of class `"mrca_age_distribution"` with distribution grid,
#'   summary statistics, requested probabilities, and the original result.
#' @export
mrca_age_distribution <- function(birth,
                                  death,
                                  time,
                                  time_step = NULL,
                                  probs = c(0.25, 0.50, 0.75),
                                  n_internal = NULL) {
  if (!is.numeric(time) || length(time) != 1L || !is.finite(time) || time < 0) {
    stop("'time' must be a finite non-negative number.", call. = FALSE)
  }
  if (is.null(time_step)) {
    time_step <- if (time > 0) time / 1000 else 1
  }
  if (!is.numeric(time_step) || length(time_step) != 1L ||
      !is.finite(time_step) || time_step <= 0) {
    stop("'time_step' must be a finite strictly positive number.", call. = FALSE)
  }
  if (!is.null(n_internal)) {
    if (!is.numeric(n_internal) || length(n_internal) != 1L || n_internal < 2) {
      stop("'n_internal' must be NULL or an integer >= 2.", call. = FALSE)
    }
    time_step <- if (time > 0) time / as.integer(n_internal) else 1
  }

  birth_fun <- .as_rate_function(birth, "birth")
  death_fun <- .as_rate_function(death, "death")
  env <- .bd_env_mrca()

  res <- env$one_mrca_time(
    T = time,
    birth_fun = birth_fun,
    death_fun = death_fun,
    dt_internal = time_step,
    probs = probs
  )

  distribution <- res$grid
  names(distribution)[names(distribution) == "s"] <- "mrca_age"

  summary <- as.data.frame(as.list(res$summary), stringsAsFactors = FALSE)
  names(summary)[names(summary) == "mean"] <- "expectation"

  out <- list(
    time = time,
    distribution = distribution,
    summary = summary,
    probs = probs,
    original_result = res
  )
  class(out) <- c("mrca_age_distribution", "list")
  out
}

#' MRCA Age Distribution Dynamics
#'
#' Computes MRCA age distribution summaries over a time grid under a birth-death
#' process.
#'
#' @inheritParams mrca_age_distribution
#' @param time_start,time_end,time_step Time grid definition.
#' @param n_internal_per_step Number of internal integration substeps per output
#'   time step.
#' @param include_zero Logical; whether to include time zero in the original
#'   dynamics computation.
#'
#' @return An object of class `"mrca_age_distribution_dynamics"` containing the
#'   summary data frame, distribution grids, probabilities, and original call.
#' @export
mrca_age_distribution_dynamics <- function(birth,
                                           death,
                                           time_start = 0,
                                           time_end,
                                           time_step,
                                           probs = c(0.25, 0.50, 0.75),
                                           n_internal_per_step = 10L,
                                           include_zero = TRUE) {
  grid <- .make_public_grid(time_start, time_end, time_step)
  if (time_end <= 0) {
    stop("'time_end' must be positive for MRCA dynamics.", call. = FALSE)
  }

  birth_fun <- .as_rate_function(birth, "birth")
  death_fun <- .as_rate_function(death, "death")
  env <- .bd_env_mrca()

  raw <- env$mrca_dynamics(
    birth = birth_fun,
    death = death_fun,
    t_max = time_end,
    dt = time_step,
    probs = probs,
    n_internal_per_step = n_internal_per_step,
    include_zero = include_zero || time_start == 0
  )

  idx <- .match_grid_indices(raw$summary$time, grid)
  raw$summary <- raw$summary[idx, , drop = FALSE]
  raw$grids <- raw$grids[idx]
  raw$call <- match.call()

  names(raw$summary)[names(raw$summary) == "mean"] <- "expectation"
  class(raw) <- c("mrca_age_distribution_dynamics", "list")
  raw
}

#' Plot MRCA Age Distribution Dynamics
#'
#' Plots MRCA age summaries over time, including median, optional interquartile
#' band, and optional expectation.
#'
#' @param x An `"mrca_age_distribution_dynamics"` object.
#' @param show_iqr Logical; whether to draw the interquartile band.
#' @param show_expectation Logical; whether to draw the expectation curve.
#' @param main,xlab,ylab Base graphics labels.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @export
plot_mrca_age_distribution_dynamics <- function(x,
                                                show_iqr = TRUE,
                                                show_expectation = TRUE,
                                                main = NA,
                                                xlab = "Time",
                                                ylab = "MRCA age",
                                                ...) {
  if (!inherits(x, "mrca_age_distribution_dynamics")) {
    stop("'x' must be the output of mrca_age_distribution_dynamics().", call. = FALSE)
  }

  d <- x$summary
  qnames <- paste0("q", x$probs * 100)

  ycols <- intersect(c("expectation", qnames), names(d))
  ylim <- range(d[, ycols, drop = FALSE], na.rm = TRUE)
  if (!all(is.finite(ylim))) {
    ylim <- c(0, max(d$time, na.rm = TRUE))
  }

  graphics::plot(d$time, d$expectation, type = "n", ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main, ...)

  q25 <- qnames[1L]
  q50 <- qnames[which.min(abs(x$probs - 0.5))]
  q75 <- qnames[length(qnames)]

  if (show_iqr && all(c(q25, q75) %in% names(d))) {
    graphics::polygon(c(d$time, rev(d$time)),
                      c(d[[q25]], rev(d[[q75]])),
                      col = grDevices::adjustcolor("grey", alpha.f = 0.35),
                      border = NA)
  }
  if (q50 %in% names(d)) {
    graphics::lines(d$time, d[[q50]], lwd = 2)
  }
  if (show_expectation && "expectation" %in% names(d)) {
    graphics::lines(d$time, d$expectation, lwd = 2, lty = 2)
  }

  graphics::legend("topleft",
                   legend = c("median", "expectation"),
                   lty = c(1, 2), lwd = 2, bty = "n")
  invisible(x)
}
