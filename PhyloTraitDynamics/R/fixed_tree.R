## Public API for fixed phylogenetic trees.

simulate_fixed_tree_brownian_realization <- function(tree,
                                                     sigma2 = 1,
                                                     time_step,
                                                     seed = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.", call. = FALSE)
  }
  if (!is.numeric(sigma2) || length(sigma2) != 1L || !is.finite(sigma2) || sigma2 < 0) {
    stop("'sigma2' must be a finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(time_step) || length(time_step) != 1L ||
      !is.finite(time_step) || time_step <= 0) {
    stop("'time_step' must be a finite strictly positive number.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  env <- .bd_env_fixed()
  df <- env$simulate_bm_on_tree(tree = tree, dt = time_step, sigma = sqrt(sigma2))
  out <- list(
    tree = tree,
    data = df,
    time = df$time,
    empirical_mean = df$empirical_mean,
    empirical_variance = df$empirical_variance,
    sigma2 = sigma2,
    time_step = time_step,
    seed = seed
  )
  class(out) <- c("fixed_tree_brownian_realization", "list")
  out
}

plot_fixed_tree_brownian_realization <- function(x,
                                                 n_bands = 10,
                                                 time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05)) {
  if (!inherits(x, "fixed_tree_brownian_realization")) {
    stop("'x' must be the output of simulate_fixed_tree_brownian_realization().", call. = FALSE)
  }
  env <- .bd_env_fixed()
  env$plot_simulation_with_tree(
    df = x$data,
    tree = x$tree,
    n_bands = n_bands,
    time_band_col = time_band_col
  )
  invisible(x)
}

fixed_tree_theoretical_summary <- function(tree,
                                           sigma2 = 1,
                                           time_start = 0,
                                           time_end = NULL,
                                           time_step) {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.", call. = FALSE)
  }
  if (!is.numeric(sigma2) || length(sigma2) != 1L || !is.finite(sigma2) || sigma2 < 0) {
    stop("'sigma2' must be a finite non-negative number.", call. = FALSE)
  }
  env <- .bd_env_fixed()
  if (is.null(time_end)) {
    time_end <- env$get_max_time(tree)
  }
  times <- .make_public_grid(time_start, time_end, time_step)

  raw <- env$compute_emp_mean_var_timeseries_general(
    tree = tree,
    times = times,
    sigma2 = sigma2
  )
  n_lineages <- .active_lineage_counts_fixed_tree(tree, times, env)

  out <- data.frame(
    time = raw$time,
    n_lineages = n_lineages,
    empirical_mean_variance = raw$var_mean,
    empirical_variance_expectation = raw$mean,
    empirical_variance_q25 = raw$q25,
    empirical_variance_q50 = raw$q50,
    empirical_variance_q75 = raw$q75,
    empirical_variance_variance = raw$var
  )
  attr(out, "tree") <- tree
  attr(out, "sigma2") <- sigma2
  attr(out, "original_result") <- raw
  out
}

plot_fixed_tree_theoretical_summary <- function(x,
                                                tree = NULL,
                                                show_tree = TRUE,
                                                main = NULL,
                                                time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05),
                                                band_col = grDevices::adjustcolor("red", alpha.f = 0.18),
                                                expectation_col = "blue",
                                                median_col = "red",
                                                mean_variance_col = "darkgreen") {
  if (!is.data.frame(x) || !all(c(
    "time", "empirical_mean_variance", "empirical_variance_expectation",
    "empirical_variance_q25", "empirical_variance_q50", "empirical_variance_q75"
  ) %in% names(x))) {
    stop("'x' must be the output of fixed_tree_theoretical_summary().", call. = FALSE)
  }
  if (is.null(tree)) {
    tree <- attr(x, "tree", exact = TRUE)
  }
  if (isTRUE(show_tree) && !inherits(tree, "phylo")) {
    stop("'tree' must be supplied when show_tree = TRUE.", call. = FALSE)
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  n_panels <- if (isTRUE(show_tree)) 3L else 2L
  if (n_panels == 3L) {
    graphics::layout(matrix(seq_len(3L), ncol = 1L), heights = c(1.2, 1, 1.2))
  } else {
    graphics::layout(matrix(seq_len(2L), ncol = 1L), heights = c(1, 1.2))
  }

  add_bands <- function() {
    u <- graphics::par("usr")
    tmin <- min(x$time, na.rm = TRUE)
    tmax <- max(x$time, na.rm = TRUE)
    width <- (tmax - tmin) / 10
    if (!is.finite(width) || width <= 0) return(invisible(NULL))
    starts <- seq(tmin, tmax, by = width)
    for (s in starts[seq_along(starts) %% 2 == 1]) {
      graphics::rect(s, u[3], min(s + width, tmax), u[4],
                     col = time_band_col, border = NA)
    }
    invisible(NULL)
  }

  if (isTRUE(show_tree)) {
    graphics::par(mar = c(0, 4, 2, 1))
    ape::plot.phylo(tree, direction = "rightwards", show.tip.label = FALSE,
                    main = main)
  }

  graphics::par(mar = c(0, 4, 1, 1))
  ymax1 <- max(x$empirical_mean_variance, na.rm = TRUE)
  if (!is.finite(ymax1) || ymax1 <= 0) ymax1 <- 1
  graphics::plot(x$time, x$empirical_mean_variance, type = "n",
                 ylim = c(0, ymax1), xlab = "", ylab = "Variance of\nempirical mean",
                 xaxt = "n")
  add_bands()
  graphics::lines(x$time, x$empirical_mean_variance,
                  col = mean_variance_col, lwd = 2)

  graphics::par(mar = c(4, 4, 1, 1))
  ymax2 <- max(x$empirical_variance_q75, x$empirical_variance_expectation,
               na.rm = TRUE)
  if (!is.finite(ymax2) || ymax2 <= 0) ymax2 <- 1
  graphics::plot(x$time, x$empirical_variance_q50, type = "n",
                 ylim = c(0, ymax2), xlab = "Time",
                 ylab = "Empirical\nvariance")
  add_bands()
  graphics::polygon(c(x$time, rev(x$time)),
                    c(x$empirical_variance_q25, rev(x$empirical_variance_q75)),
                    col = band_col, border = NA)
  graphics::lines(x$time, x$empirical_variance_q25, col = median_col, lty = 2)
  graphics::lines(x$time, x$empirical_variance_q75, col = median_col, lty = 2)
  graphics::lines(x$time, x$empirical_variance_q50, col = median_col, lwd = 2)
  graphics::lines(x$time, x$empirical_variance_expectation,
                  col = expectation_col, lwd = 2)
  graphics::legend("topleft",
                   legend = c("Expectation", "Median", "Quartiles"),
                   col = c(expectation_col, median_col, median_col),
                   lty = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n")

  invisible(x)
}
