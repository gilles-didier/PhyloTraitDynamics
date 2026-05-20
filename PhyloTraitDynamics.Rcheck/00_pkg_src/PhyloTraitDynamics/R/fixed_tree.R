# Brownian empirical variance on a fixed phylogenetic tree.
#
# The public API is restricted to:
#   - fixed_tree_compute_theoretical_summary()
#   - fixed_tree_plot_theoretical_summary()
#   - fixed_tree_simulate_brownian_realization()
#   - fixed_tree_plot_brownian_realization()
#
# The remaining functions are internal helpers adapted from the original script.

#' Theoretical Summary on a Fixed Tree
#'
#' Computes theoretical summaries of empirical mean variance and empirical
#' variance for Brownian traits on a fixed phylogenetic tree.
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param sigma2 Brownian variance parameter.
#' @param time_end End of the time grid. If `NULL`, the maximum tree time is
#'   used.
#' @param time_step Time step used to define the time grid.
#'
#' @return A data frame with time, lineage counts, empirical mean variance,
#'   empirical variance expectation, empirical variance quartiles, and empirical
#'   variance variance. The tree, `sigma2`, and original result are attached as
#'   attributes.
#'
#' @export
fixed_tree_compute_theoretical_summary <- function(tree,
                                                   sigma2 = 1,
                                                   time_end = NULL,
                                                   time_step) {
  time_start <- 0
  if (is.null(time_end)) {
    time_end <- .fixed_tree_get_max_time(tree)
  }

  times <- seq(time_start, time_end, by = time_step)

  original_result <- .fixed_tree_compute_emp_mean_var_timeseries_general(
    tree = tree,
    times = times,
    sigma2 = sigma2
  )

  n_lineages <- .fixed_tree_compute_lineage_counts(tree = tree, times = times)

  out <- data.frame(
    time = original_result$time,
    n_lineages = n_lineages,
    empirical_mean_variance = original_result$var_mean,
    empirical_variance_expectation = original_result$mean,
    empirical_variance_q25 = original_result$q25,
    empirical_variance_q50 = original_result$q50,
    empirical_variance_q75 = original_result$q75,
    empirical_variance_variance = original_result$var,
    check.names = FALSE
  )

  attr(out, "tree") <- tree
  attr(out, "sigma2") <- sigma2
  attr(out, "original_result") <- original_result
  class(out) <- c("fixed_tree_theoretical_summary", class(out))

  out
}

#' Plot a Fixed-Tree Theoretical Summary
#'
#' Plots fixed-tree theoretical summaries for empirical mean variance and
#' empirical variance, optionally with the tree.
#'
#' @param x A `fixed_tree_theoretical_summary` object, as returned by
#'   `fixed_tree_compute_theoretical_summary()`.
#' @param tree Optional phylogenetic tree of class `"phylo"`. If `NULL`, the
#'   tree attached to `x` is used when available.
#' @param show_tree Logical; whether to show the tree panel.
#' @param main Optional main title.
#' @param time_band_col Color used for time bands.
#' @param band_col Color used for the empirical-variance interquartile band.
#' @param expectation_col Color used for the empirical-variance expectation.
#' @param median_col Color used for the empirical-variance median and quartiles.
#' @param mean_variance_col Color used for the empirical-mean variance.
#' @param cex Character expansion factor passed to base graphics parameters.
#'
#' @return Invisibly returns `x`.
#'
#' @export
fixed_tree_plot_theoretical_summary <- function(
  x,
  tree = NULL,
  show_tree = TRUE,
  main = NULL,
  time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05),
  band_col = grDevices::adjustcolor("red", alpha.f = 0.18),
  expectation_col = "blue",
  median_col = "red",
  mean_variance_col = "darkgreen",
  cex = 1
) {
  theory <- .fixed_tree_extract_original_theory(x)

  if (is.null(tree)) {
    tree <- attr(x, "tree")
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  bands <- .fixed_tree_compute_band_limits(0, max(theory$time), n_bands = 10)

  if (isTRUE(show_tree) && !is.null(tree)) {
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("Package 'phytools' is required when show_tree = TRUE.")
    }

    graphics::layout(matrix(c(1, 2, 3), ncol = 1),
                     heights = c(1, 1.25, 2))
	graphics::par(cex = cex)
    graphics::par(mar = c(0, 4, 1, 0.1))
    phytools::plotTree(tree,
                       direction = "rightwards",
                       ftype = "off",
                       x.lim = c(0, max(theory$time)),
                       mar = c(0, 4, 1, 0.1),
                       lwd = 1)
    .fixed_tree_add_time_bands(bands, col = time_band_col)
    if (!is.null(main)) {
      graphics::title(main = main)
    }

    graphics::par(mar = c(0, 4, 0, 0.1))
    .fixed_tree_plot_empmean_var(theory, mean_variance_col = mean_variance_col)
    .fixed_tree_add_time_bands(bands, col = time_band_col)

    graphics::par(mar = c(4, 4, 0, 0.1))
    .fixed_tree_plot_empvar_quartiles(
      theory,
      band_col = band_col,
      expectation_col = expectation_col,
      median_col = median_col
    )
    .fixed_tree_add_time_bands(bands, col = time_band_col)
  } else {
    graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(1.25, 2))

    graphics::par(mar = c(0, 4, 1, 0.1))
    .fixed_tree_plot_empmean_var(theory, mean_variance_col = mean_variance_col)
    .fixed_tree_add_time_bands(bands, col = time_band_col)
    if (!is.null(main)) {
      graphics::title(main = main)
    }

    graphics::par(mar = c(4, 4, 0, 0.1))
    .fixed_tree_plot_empvar_quartiles(
      theory,
      band_col = band_col,
      expectation_col = expectation_col,
      median_col = median_col
    )
    .fixed_tree_add_time_bands(bands, col = time_band_col)
  }

  invisible(x)
}

#' Simulate Brownian Traits on a Fixed Tree
#'
#' Simulates one Brownian trait realization on a fixed phylogenetic tree and
#' records empirical mean and empirical variance through time.
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param sigma2 Brownian variance parameter.
#' @param time_step Time step used for the simulation output.
#' @param seed Optional random seed.
#'
#' @return An object of class `"fixed_tree_brownian_realization"` containing
#'   the input tree, simulation data frame, time series, parameters, and seed.
#'
#' @export
fixed_tree_simulate_brownian_realization <- function(tree,
                                                     sigma2 = 1,
                                                     time_step,
                                                     seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  df <- .fixed_tree_simulate_bm_on_tree(
    tree = tree,
    dt = time_step,
    sigma = sqrt(sigma2)
  )

  structure(
    list(
      tree = tree,
      data = df,
      time = df$time,
      sigma2 = sigma2,
      time_step = time_step,
      parameters = list(sigma2 = sigma2, time_step = time_step),
      seed = seed
    ),
    class = "fixed_tree_brownian_realization"
  )
}

#' Plot a Fixed-Tree Brownian Realization
#'
#' Plots the simulated Brownian realization together with the fixed tree using
#' the original plotting routine.
#'
#' @param x A `fixed_tree_brownian_realization` object.
#' @param n_bands Number of time bands shown by the original plotting routine.
#' @param time_band_col Color used for time bands.
#' @param cex Character expansion factor passed to base graphics parameters.
#'
#' @return Invisibly returns `x`.
#'
#' @export
fixed_tree_plot_brownian_realization <- function(
  x,
  n_bands = 10,
  time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05),
  cex = 1
) {
  if (!inherits(x, "fixed_tree_brownian_realization")) {
    stop("'x' must be a 'fixed_tree_brownian_realization' object.")
  }

    df <- x$data
    tree <- x$tree
	branch_cols <- grep("^branch_", colnames(df))
	nedge <- nrow(tree$edge)

	branch_color <- c("blue", "red", "green", "orange", "cyan", "brown", "darkgreen", "deeppink")
	branch_color <- rep(branch_color, length.out = nedge)

	bands <- .fixed_tree_compute_band_limits(0, max(df$time), n_bands)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))

  graphics::layout(matrix(1:4, ncol = 1), heights = c(2, 3, 2, 3))
    graphics::par(cex = cex)
	## Panel 1: tree
	graphics::par(mar = c(0, 4, 1, 0.1))

	graphics::plot(tree,
	     direction = "rightwards",
	     show.tip.label = FALSE,
	     edge.color = branch_color,
	     edge.width = 3)

	.fixed_tree_add_time_bands(bands, col = time_band_col)

	## Panel 2: branch values
	graphics::par(mar = c(0, 4, 0, 0.1))

	graphics::plot(df$time, df$branch_1,
	     type = "n",
	     ylim = range(df[, branch_cols]*2, na.rm = TRUE),
	     xlab = "",
	     ylab = "",
	     bty = "n",
	     xaxt = "n",
	     yaxt = "n",
	     mgp = c(1, 0.6, 0))

	.fixed_tree_add_time_bands(bands, col = time_band_col)

	graphics::mtext("Brownian\ntrait", side = 2, line = 0, cex = cex)
	graphics::lines(df$time, df$branch_1,
	      col = grDevices::adjustcolor(branch_color[1], alpha.f = 0.7),
	      lwd = 1)

	if (nedge >= 2) {
		for (j in 2:nedge) {
			graphics::lines(df$time, df[[paste0("branch_", j)]],
			      col = grDevices::adjustcolor(branch_color[j], alpha.f = 0.7),
			      lwd = 1)
		}
	}

  ## Panel 3: empirical mean
  graphics::par(mar = c(0, 4, 0, 0.1))

  graphics::plot(df$time, df[, "empirical_mean"],
       type = "n",
       ylim = range(df[, "empirical_mean"]*2, na.rm = TRUE),
       xlab = "",
       ylab = "",
       bty = "n",
       xaxt = "n",
       yaxt = "n",
       mgp = c(1, 0.6, 0))

  .fixed_tree_add_time_bands(bands, col = time_band_col)

  graphics::mtext("Empirical\nmean", side = 2, line = 0, cex = cex)
  graphics::lines(df$time, df[, "empirical_mean"],
        col = "black",
        lwd = 1)

	## Panel 4: empirical variance
	graphics::par(mar = c(4, 4, 0, 0.1))

	graphics::plot(df$time, df[, "empirical_variance"],
	     type = "n",
	     ylim = range(df[, "empirical_variance"]*2, na.rm = TRUE),
	     xlab = "Time",
	     ylab = "",
	     bty = "n",
	     xaxt = "s",
	     yaxt = "n",
	     mgp = c(1, 0.6, 0))

	.fixed_tree_add_time_bands(bands, col = time_band_col)
	graphics::mtext("Empirical\nvariance", side = 2, line = 0, cex = cex)
	graphics::lines(df$time, df[, "empirical_variance"],
	      col = "black",
	      lwd = 1)

  invisible(x)
}

.fixed_tree_compute_node_depths <- function(tree) {
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  n_tip <- length(tree$tip.label)
  depths <- numeric(n_tip + tree$Nnode)
  root <- n_tip + 1
  depths[root] <- 0

  for (intern in seq.int(root, root + tree$Nnode - 1)) {
    idx <- which(tree$edge[, 1] == intern)
    if (length(idx) == 0) next
    for (e in idx) {
      child <- tree$edge[e, 2]
      depths[child] <- depths[intern] + tree$edge.length[e]
    }
  }
  depths
}

.fixed_tree_get_max_time <- function(tree) {
  depths <- .fixed_tree_compute_node_depths(tree)
  n_tip <- length(tree$tip.label)
  max(depths[seq_len(n_tip)])
}

.fixed_tree_get_mrca <- function(tree, a, b) {
  if (a != b) {
    n_tip <- length(tree$tip.label)
    if (a <= n_tip) {
      a <- tree$edge[tree$edge[, 2] == a, 1]
    }
    if (b <= n_tip) {
      b <- tree$edge[tree$edge[, 2] == b, 1]
    }
    while (a != b) {
      if (a < b) {
        b <- tree$edge[tree$edge[, 2] == b, 1]
      } else {
        a <- tree$edge[tree$edge[, 2] == a, 1]
      }
    }
  }
  a
}

.fixed_tree_compute_empvar_quartiles <- function(C) {
  if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("Package 'CompQuadForm' is required for theoretical quartiles.")
  }

  n <- nrow(C)
  P <- diag(n) - matrix(1 / n, n, n)
  H <- chol(C)
  evals <- eigen(H %*% P %*% t(H), symmetric = TRUE, only.values = TRUE)$values
  evals <- evals[evals > 1e-12]

  if (length(evals) == 0) {
    return(c(q25 = 0, q50 = 0, q75 = 0, mean = 0, var = 0))
  }

  upper <- sum(evals) * max(stats::qchisq(0.9999, df = length(evals)), 10)

  f25 <- function(x) CompQuadForm::davies(x, lambda = evals)$Qq - 0.25
  f50 <- function(x) CompQuadForm::davies(x, lambda = evals)$Qq - 0.50
  f75 <- function(x) CompQuadForm::davies(x, lambda = evals)$Qq - 0.75

  q25 <- stats::uniroot(f25, lower = 0, upper = upper)$root / (n - 1)
  q50 <- stats::uniroot(f50, lower = 0, upper = upper)$root / (n - 1)
  q75 <- stats::uniroot(f75, lower = 0, upper = upper)$root / (n - 1)
  mean <- sum(evals) / (n - 1)
  var  <- 2 * sum(evals^2) / (n - 1)^2

  c(q25 = q25, q50 = q50, q75 = q75, mean = mean, var = var)
}

.fixed_tree_compute_emp_mean_var_timeseries_general <- function(
  tree,
  times = seq(0, .fixed_tree_get_max_time(tree), by = 0.1),
  sigma2 = 1
) {
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  depths <- .fixed_tree_compute_node_depths(tree)
  total <- max(tree$edge)

  CTotal <- matrix(0, total, total)
  for (i in seq_len(total)) {
    for (j in i:total) {
      cij <- depths[.fixed_tree_get_mrca(tree, i, j)]
      CTotal[i, j] <- cij
      CTotal[j, i] <- cij
    }
  }

  edge_start <- depths[tree$edge[, 1]]
  edge_end   <- depths[tree$edge[, 2]]

  results <- vector("list", length(times))

  for (k in seq_along(times)) {
    t <- times[k]

    active_edges <- which(edge_start < t & edge_end >= t)
    n_now <- length(active_edges)

    if (n_now < 2) {
      results[[k]] <- data.frame(time = t, q25 = 0, q50 = 0, q75 = 0, mean = 0, var = 0, var_mean = 0)
      next
    }

    active_nodes <- tree$edge[active_edges, 2]
    C <- matrix(0, n_now, n_now)

    for (i in seq_len(n_now)) {
      C[i, i] <- t
      if (i < n_now) {
        for (j in (i + 1):n_now) {
          C[i, j] <- CTotal[active_nodes[i], active_nodes[j]]
          C[j, i] <- C[i, j]
        }
      }
    }

    qu <- .fixed_tree_compute_empvar_quartiles(sigma2 * C)
    results[[k]] <- data.frame(
      time = t,
      q25  = unname(qu["q25"]),
      q50  = unname(qu["q50"]),
      q75  = unname(qu["q75"]),
      mean = unname(qu["mean"]),
      var  = unname(qu["var"]),
      var_mean = sum(C[ , ])/(n_now**2)
    )
  }

  do.call(rbind, results)
}

.fixed_tree_compute_lineage_counts <- function(tree, times) {
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  depths <- .fixed_tree_compute_node_depths(tree)

  edge_start <- depths[tree$edge[, 1]]
  edge_end   <- depths[tree$edge[, 2]]

  vapply(times, function(t) {
    length(which(edge_start < t & edge_end >= t))
  }, integer(1))
}

.fixed_tree_compute_band_limits <- function(tmin, tmax, n_bands = 10) {
  width_raw <- (tmax - tmin) / n_bands
  step <- 10 ^ floor(log10(width_raw))
  if(step == width_raw) {
	  width = width_raw
  } else {
	  width = step*10
	  min = abs(width_raw-width)
	  if(abs(width_raw-width*0.5)<min) {
		  width = width*0.5
		  min = abs(width_raw-width)
	 }
	  if(abs(width_raw-width*0.5)<min) {
		  width = width*0.5
		  min = abs(width_raw-width)
	 }
	  if(abs(width_raw-step)<min) {
		  width = step
		  min = abs(width_raw-width)
	 }
  }
  first_start <- floor(tmin / width) * width
  starts <- seq(first_start, tmax, by = width)
  ends <- starts + width / 2
  ends[length(ends)] <- min(ends[length(ends)], tmax)
  data.frame(start = starts, end = ends)
}

.fixed_tree_add_time_bands <- function(bands, col = grDevices::adjustcolor("grey", alpha.f = 0.2)) {
  usr <- graphics::par("usr")
  for (i in seq_len(nrow(bands))) {
    graphics::rect(bands$start[i], usr[3], bands$end[i], usr[4], col = col, border = NA)
  }
}

.fixed_tree_extract_original_theory <- function(x) {
  original_result <- attr(x, "original_result")
  if (!is.null(original_result)) {
    return(original_result)
  }

  if (all(c("time", "q25", "q50", "q75", "mean", "var", "var_mean") %in% colnames(x))) {
    return(x)
  }

  data.frame(
    time = x$time,
    q25 = x$empirical_variance_q25,
    q50 = x$empirical_variance_q50,
    q75 = x$empirical_variance_q75,
    mean = x$empirical_variance_expectation,
    var = x$empirical_variance_variance,
    var_mean = x$empirical_mean_variance,
    check.names = FALSE
  )
}

.fixed_tree_plot_empvar_quartiles <- function(df,
                                              band_col = grDevices::rgb(1, 0, 0, alpha = 0.2),
                                              expectation_col = "blue",
                                              median_col = "red") {
	if(!all(c("time", "q25", "q50", "q75", "mean", "var") %in% colnames(df))) {
		stop("The dataframe must contain columns : time, q25, q50, q75, mean")
	}
	# Median
	graphics::plot(df$time, df$q50, type = "l", col = median_col, lty = 1, lwd = 1, ylim = range(c(df$q25, df$q75, df$mean, df$q50)), axes = FALSE, xlab = "Time", ylab = "Empirical variance\ndistribution", mgp = c(2, 0.6, 0))
    graphics::Axis(side=1)
    graphics::Axis(side=2)

	# Quartiles
	graphics::lines(df$time, df$q25, col = median_col, lty = 2, lwd = 1)
	graphics::lines(df$time, df$q75, col = median_col, lty = 2, lwd = 1)

	# Shade between q25 et q75
	graphics::polygon(c(df$time, rev(df$time)),
	c(df$q25, rev(df$q75)),
	col = band_col, border = NA)

	# Mean
	graphics::lines(df$time, df$mean, col = expectation_col, lty = 1, lwd = 0.5)

	graphics::legend("topleft", legend = c("Expectation", "Median", "Quartiles"),
	col = c(expectation_col, median_col, median_col), lty = c(1,1,2), lwd = c(0.5, 1.5, 1),
	bty = "n")
}

.fixed_tree_plot_empmean_var <- function(df, mean_variance_col = "darkgreen") {
	if(!all(c("time", "q25", "q50", "q75", "mean", "var") %in% colnames(df))) {
		stop("The dataframe must contain columns : time, q25, q50, q75, mean")
	}
	# Median
	graphics::plot(df$time, df$var_mean, type = "l", col = mean_variance_col, lty = 1, lwd = 1.5, ylim = range(c(df$var_mean)), axes = FALSE, xlab = NULL, ylab = "Empirical mean\ndistribution", mgp = c(2, 0.6, 0))
    graphics::Axis(side=2)

	graphics::legend("topleft", legend = c("Variance"),
	col = c(mean_variance_col), lty = c(1), lwd = c(1.5),
	bty = "n")
}

.fixed_tree_simulate_bm_on_tree <- function(tree, dt = 0.01, sigma = 1) {
	depth <- .fixed_tree_compute_node_depths(tree)
	Tmax <- max(depth)
	time <- seq(0, Tmax, by = dt)
	nt <- length(time)
	root <- length(tree$tip.label)+1
	value <- numeric(length(depth))
	value[root] <- 0.
	mat <- matrix(NA, nt, nrow(tree$edge))
	colnames(mat) <- paste0("branch_", seq_len(nrow(tree$edge)))
	for(intern in seq(root,root+tree$Nnode-1)) {
		t_start <- depth[intern]
		for(edge in which(tree$edge[,1] == intern)) {
			t_end <- depth[tree$edge[edge, 2]]
			idx <- which(time >= t_start & time <= t_end)
			n <- length(idx)
			increments <- stats::rnorm(n - 1, mean = 0, sd = sigma * sqrt(dt))
			path <- cumsum(c(value[intern], increments))
			mat[idx, edge] <- path
			value[tree$edge[edge, 2]] <- path[length(path)]
		}
	}
	df <- data.frame(time = time, mat)
	
	branch_values <- df[, grep("^branch_", colnames(df)), drop = FALSE]

	df$empirical_variance <- apply(branch_values, 1, function(x) {
	  x <- x[!is.na(x)]
	  if (length(x) < 2L) return(0)
	  stats::var(x)
	})

	df$empirical_mean <- apply(branch_values, 1, function(x) {
	  x <- x[!is.na(x)]
	  if (length(x) == 0L) return(NA_real_)
	  mean(x)
	})
  	df
}
