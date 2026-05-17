# Brownian empirical variance on a phylogenetic tree
# Generic version: works for trees with or without extinctions.
#
# Main functions:
#   - simulate_empvar_bm(tree, n, sigma2, step)
#   - compute_empvar_timeseries_general(tree, times, sigma2)
#   - plot_empvar_bm(tree, sim)
#
# Here sigma2 is the Brownian variance rate:
#   Var(B(t+h)-B(t)) = sigma2 * h

compute_node_depths <- function(tree) {
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

get_max_time <- function(tree) {
  depths <- compute_node_depths(tree)
  n_tip <- length(tree$tip.label)
  max(depths[seq_len(n_tip)])
}

getMRCA_OK <- function(tree, a, b) {
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

compute_empvar_quartiles <- function(C) {
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

compute_empvar_timeseries_general <- function(tree,
                                              times = seq(0, get_max_time(tree), by = 0.1),
                                              sigma2 = 1) {
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  depths <- compute_node_depths(tree)
  total <- max(tree$edge)

  CTotal <- matrix(0, total, total)
  for (i in seq_len(total)) {
    for (j in i:total) {
      cij <- depths[getMRCA_OK(tree, i, j)]
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
      results[[k]] <- data.frame(time = t, q25 = 0, q50 = 0, q75 = 0, mean = 0, var = 0)
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

    qu <- compute_empvar_quartiles(sigma2 * C)
    results[[k]] <- data.frame(
      time = t,
      q25  = unname(qu["q25"]),
      q50  = unname(qu["q50"]),
      q75  = unname(qu["q75"]),
      mean = unname(qu["mean"]),
      var  = unname(qu["var"])
    )
  }

  do.call(rbind, results)
}

compute_emp_mean_var_timeseries_general <- function(tree,
                                              times = seq(0, get_max_time(tree), by = 0.1),
                                              sigma2 = 1) {
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  depths <- compute_node_depths(tree)
  total <- max(tree$edge)

  CTotal <- matrix(0, total, total)
  for (i in seq_len(total)) {
    for (j in i:total) {
      cij <- depths[getMRCA_OK(tree, i, j)]
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

    qu <- compute_empvar_quartiles(sigma2 * C)
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

simulate_empvar_bm <- function(tree, n, sigma2, step = 0.01, seed = NULL, theory = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.")
  }
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("'n' must be a positive integer.")
  }
  if (!is.numeric(sigma2) || length(sigma2) != 1 || sigma2 < 0) {
    stop("'sigma2' must be a nonnegative number.")
  }
  if (!is.numeric(step) || length(step) != 1 || step <= 0) {
    stop("'step' must be a positive number.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  tree <- ape::reorder.phylo(tree, order = "cladewise")
  depths <- compute_node_depths(tree)
  Tmax <- get_max_time(tree)

  times <- seq(0, Tmax, by = step)
  if (tail(times, 1) < Tmax) {
    times <- c(times, Tmax)
  }

  n_tip <- length(tree$tip.label)
  root <- n_tip + 1
  n_edge <- nrow(tree$edge)
  n_total <- n_tip + tree$Nnode

  edge_start <- depths[tree$edge[, 1]]
  edge_end   <- depths[tree$edge[, 2]]
  active_edges <- lapply(times, function(t) which(edge_start < t & edge_end >= t))

  simulate_once <- function() {
    node_value <- numeric(n_total)
    node_value[root] <- 0
    edge_value <- matrix(NA_real_, nrow = length(times), ncol = n_edge)

    for (intern in seq.int(root, root + tree$Nnode - 1)) {
      child_edges <- which(tree$edge[, 1] == intern)
      if (length(child_edges) == 0) next

      for (e in child_edges) {
        t0 <- depths[tree$edge[e, 1]]
        t1 <- depths[tree$edge[e, 2]]

        local_times <- c(t0, times[times > t0 & times < t1], t1)
        if (length(local_times) == 1) {
          x <- node_value[intern]
        } else {
          inc <- stats::rnorm(length(local_times) - 1,
                              mean = 0,
                              sd = sqrt(sigma2 * diff(local_times)))
          x <- node_value[intern] + c(0, cumsum(inc))
        }

        idx_global <- which(times > t0 & times <= t1)
        if (length(idx_global) > 0) {
          pos <- match(times[idx_global], local_times)
          edge_value[idx_global, e] <- x[pos]
        }

        node_value[tree$edge[e, 2]] <- x[length(x)]
      }
    }

    empvar <- numeric(length(times))
    for (k in seq_along(times)) {
      ee <- active_edges[[k]]
      if (length(ee) < 2) {
        empvar[k] <- 0
      } else {
        empvar[k] <- stats::var(edge_value[k, ee])
      }
    }
    empvar
  }

  sim_mat <- replicate(as.integer(n), simulate_once())
  if (is.null(dim(sim_mat))) {
    sim_mat <- matrix(sim_mat, ncol = 1)
  }
  colnames(sim_mat) <- paste0("sim_", seq_len(ncol(sim_mat)))

  sim_df <- data.frame(time = times, sim_mat, check.names = FALSE)

  theo_df <- NULL
  if (isTRUE(theory)) {
    theo_df <- compute_emp_mean_var_timeseries_general(tree, times = times, sigma2 = sigma2)
  }

  structure(
    list(
      time = times,
      simulations = sim_mat,
      sim_df = sim_df,
      theory = theo_df,
      sigma2 = sigma2,
      step = step,
      tree = tree
    ),
    class = "empvar_bm_sim"
  )
}

compute_band_limits <- function(tmin, tmax, n_bands = 10) {
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

add_time_bands <- function(bands, col = grDevices::adjustcolor("grey", alpha.f = 0.2)) {
  usr <- par("usr")
  for (i in seq_len(nrow(bands))) {
    rect(bands$start[i], usr[3], bands$end[i], usr[4], col = col, border = NA)
  }
}

plot_empvar_quartiles <- function(df) {
	if(!all(c("time", "q25", "q50", "q75", "mean", "var") %in% colnames(df))) {
		stop("The dataframe must contain columns : time, q25, q50, q75, mean")
	}
	cexA <- 1
	par(cex.axis = cexA*0.9, cex.lab = cexA, cex.main = cexA)
	# Median
	plot(df$time, df$q50, type = "l", col = "red", lty = 1, lwd = 1, ylim = range(c(df$q25, df$q75, df$mean, df$q50)), axes = FALSE, xlab = "Time", ylab = "Distribution of\n the empirical variance", cex.lab = cexA*1.1, cex.axis = cexA*1.1, mgp = c(2, 0.6, 0))
Axis(side=1, cex.axis = cexA*1.1)
Axis(side=2, cex.axis = cexA*1.1)

	# Quartiles
	lines(df$time, df$q25, col = "red", lty = 2, lwd = 1)
	lines(df$time, df$q75, col = "red", lty = 2, lwd = 1)

	# Shade between q25 et q75
	polygon(c(df$time, rev(df$time)),
	c(df$q25, rev(df$q75)),
	col = rgb(1,0,0,alpha = 0.2), border = NA)

	# Mean
	lines(df$time, df$mean, col = "blue", lty = 1, lwd = 0.5)

	legend("topleft", legend = c("Expectation", "Median", "Quartiles"),
	col = c("blue", "red", "red"), lty = c(1,1,2), lwd = c(0.5, 1.5, 1),
	bty = "n", cex = cexA*1.1)
}

plot_empvar_with_tree_bm <- function(df, tree, sim, show_tree = TRUE, n_bands = 10,
                           time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05),
                           sim_col = grDevices::rgb(0, 0, 0, alpha = 0.08),
                           band_col = grDevices::rgb(1, 0, 0, alpha = 0.18),
                           median_col = "red",
                           mean_col = "blue") {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.")
  }
  if (!is.list(sim) || is.null(sim$time) || is.null(sim$simulations)) {
    stop("'sim' must be the output of simulate_empvar_bm().")
  }
  times <- sim$time
  S <- sim$simulations
  theory <- sim$theory
  ymax <- max(S, na.rm = TRUE)
  if (!is.null(theory)) {
    ymax <- max(ymax, theory$q75, theory$mean, na.rm = TRUE)
  }
  ymax <- max(ymax, 1e-12)
	bands <- compute_band_limits(0, max(sim$time), n_bands)
#	par(mfrow = c(2,1))
nf <- layout(matrix(c(1,2),ncol=1), widths=c(4,4,4), heights=c(1,2), TRUE) 
	par(mar = c(0,4,1,0.1))
	plotTree(tree, direction = "rightwards", ftype = "off", x.lim = c(0, max(sim$time)), mar = c(0,4,1,0.1), lwd = 1)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=time_band_col, border=NA)
	}
	par(mar = c(4,4,0,0.1))
	plot_empvar_quartiles(theory)
#	matlines(times, S, lty = 1, col = sim_col)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=time_band_col, border=NA)
	}
}

plot_empmean_var <- function(df) {
	if(!all(c("time", "q25", "q50", "q75", "mean", "var") %in% colnames(df))) {
		stop("The dataframe must contain columns : time, q25, q50, q75, mean")
	}
	cexA <- 1
	par(cex.axis = cexA*0.9, cex.lab = cexA, cex.main = cexA)
	# Median
	plot(df$time, df$var_mean, type = "l", col = "darkgreen", lty = 1, lwd = 1.5, ylim = range(c(df$var_mean)), axes = FALSE, xlab = NULL, ylab = "Distribution of\n the empirical mean", cex.lab = cexA*1.1, cex.axis = cexA*1.1, mgp = c(2, 0.6, 0))
Axis(side=2, cex.axis=cexA*1.1)

	legend("topleft", legend = c("Variance"),
	col = c("darkgreen"), lty = c(1), lwd = c(1.5),
	bty = "n", cex = cexA*1.1)
}

plot_emp_mean_var_with_tree_bm <- function(df, tree, sim, show_tree = TRUE, n_bands = 10,
                           time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05),
                           sim_col = grDevices::rgb(0, 0, 0, alpha = 0.08),
                           band_col = grDevices::rgb(1, 0, 0, alpha = 0.18),
                           median_col = "red",
                           mean_col = "blue") {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.")
  }
  if (!is.list(sim) || is.null(sim$time) || is.null(sim$simulations)) {
    stop("'sim' must be the output of simulate_empvar_bm().")
  }
  times <- sim$time
  S <- sim$simulations
  theory <- sim$theory
  ymax <- max(S, na.rm = TRUE)
  if (!is.null(theory)) {
    ymax <- max(ymax, theory$q75, theory$mean, na.rm = TRUE)
  }
  ymax <- max(ymax, 1e-12)
	bands <- compute_band_limits(0, max(sim$time), n_bands)
nf <- layout(matrix(c(1,2,3),ncol=1), widths=c(6,6,6), heights=c(1,1.25,2), TRUE) 
	par(mar = c(0,4,1,0.1))
	plotTree(tree, direction = "rightwards", ftype = "off", x.lim = c(0, max(sim$time)), mar = c(0,4,1,0.1), lwd = 1)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=time_band_col, border=NA)
	}
	par(mar = c(0,4,0,0.1))
	plot_empmean_var(theory)
#	matlines(times, S, lty = 1, col = sim_col)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=time_band_col, border=NA)
	}
	par(mar = c(4,4,0,0.1))
	plot_empvar_quartiles(theory)
#	matlines(times, S, lty = 1, col = sim_col)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=time_band_col, border=NA)
	}
}


plot_empvar_with_tree_bm_bis <- function(df, tree, sim, show_tree = TRUE, n_bands = 10,
                           time_band_col = grDevices::adjustcolor("grey", alpha.f = 0.2),
                           sim_col = grDevices::rgb(0, 1, 0, alpha = 0.08),
                           band_col = grDevices::rgb(1, 0, 0, alpha = 0.18),
                           median_col = "red",
                           mean_col = "blue") {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.")
  }
  if (!is.list(sim) || is.null(sim$time) || is.null(sim$simulations)) {
    stop("'sim' must be the output of simulate_empvar_bm().")
  }
  times <- sim$time
  S <- sim$simulations
  theory <- sim$theory
  ymax <- max(S, na.rm = TRUE)
  if (!is.null(theory)) {
    ymax <- max(ymax, theory$q75, theory$mean, na.rm = TRUE)
  }
  ymax <- max(ymax, 1e-12)
	bands <- compute_band_limits(0, max(sim$time), n_bands)
	par(mfrow = c(3,1))
	par(mar = c(0,4,1,0.1))
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=make.transparent("grey",0.2), border=NA)
	}
	plotTree(tree, direction = "rightwards", ftype = "off", x.lim = c(0, max(sim$time)), mar = c(0,4,1,0.1), lwd = 0.5)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=make.transparent("grey",0.2), border=NA)
	}
	par(mar = c(4,4,0,0.1))
	plot_empvar_quartiles(theory)
	matlines(times, S, lty = 1, col = sim_col)
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=make.transparent("grey",0.2), border=NA)
	}
	par(mar = c(4,4,0,0.1))
	plot(theory$time, sqrt(theory$var), type = "l", col = "blue", lty = 1, lwd = 1, ylim = range(sqrt(theory$var)), xlab = "Time", ylab = "Empirical variance", cex.lab = 1.75, cex.axis = 1.25, mgp = c(2, 0.6, 0))
	for(i in seq_len(nrow(bands))) {
		rect(bands$start[i], par()$usr[3], bands$end[i], par()$usr[4], ,col=make.transparent("grey",0.2), border=NA)
	}
}

plot_empvar_bm <- function(tree, sim, show_tree = TRUE,
                           n_bands = 10,
                           time_band_col = grDevices::adjustcolor("grey", alpha.f = 0.2),
                           sim_col = grDevices::rgb(0, 0, 0, alpha = 0.08),
                           band_col = grDevices::rgb(1, 0, 0, alpha = 0.18),
                           median_col = "red",
                           mean_col = "blue") {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be of class 'phylo'.")
  }
  if (!is.list(sim) || is.null(sim$time) || is.null(sim$simulations)) {
    stop("'sim' must be the output of simulate_empvar_bm().")
  }

  times <- sim$time
  S <- sim$simulations
  theory <- sim$theory
  bands <- compute_band_limits(0, max(times), n_bands = n_bands)

  ymax <- max(S, na.rm = TRUE)
  if (!is.null(theory)) {
    ymax <- max(ymax, theory$q75, theory$mean, na.rm = TRUE)
  }
  ymax <- max(ymax, 1e-12)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (isTRUE(show_tree)) {
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("Package 'phytools' is required when show_tree = TRUE.")
    }
    layout(matrix(c(1, 2), nrow = 2), heights = c(1, 2))

    par(mar = c(0, 4, 1, 0.5))
    phytools::plotTree(tree, direction = "rightwards", ftype = "off",
                       x.lim = c(0, max(times)), lwd = 0.5)
    add_time_bands(bands, col = time_band_col)

    par(mar = c(4, 4, 0, 0.5))
  } else {
    par(mar = c(4, 4, 1, 0.5))
  }

  plot(times, S[, 1], type = "n", xlab = "Time", ylab = "Empirical variance",
       ylim = c(0, ymax), xaxs = "i")
  add_time_bands(bands, col = time_band_col)
  matlines(times, S, lty = 1, col = sim_col)

  if (!is.null(theory)) {
    polygon(c(times, rev(times)), c(theory$q25, rev(theory$q75)),
            col = band_col, border = NA)
    lines(times, theory$q25, col = median_col, lty = 2, lwd = 1)
    lines(times, theory$q75, col = median_col, lty = 2, lwd = 1)
    lines(times, theory$q50, col = median_col, lty = 1, lwd = 1.5)
    lines(times, theory$mean, col = mean_col, lty = 1, lwd = 2)

    legend("topleft",
           legend = c("Expectation", "Median", "Quartiles", "Simulations"),
           col = c(mean_col, median_col, median_col, sim_col),
           lty = c(1, 1, 2, 1),
           lwd = c(2, 1.5, 1, 1),
           bty = "n")
  } else {
    legend("topleft", legend = "Simulations", col = sim_col, lty = 1, bty = "n")
  }

  invisible(NULL)
}

simulate_bm_on_tree <- function(tree, dt = 0.01, sigma = 1) {
	depth <- compute_node_depths(tree)
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
			increments <- rnorm(n - 1, mean = 0, sd = sigma * sqrt(dt))
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
	  var(x)
	})

	df$empirical_mean <- apply(branch_values, 1, function(x) {
	  x <- x[!is.na(x)]
	  if (length(x) == 0L) return(NA_real_)
	  mean(x)
	})
  	df
}


plot_simulation_with_tree <- function(df, tree, n_bands = 10,
                                      time_band_col = grDevices::adjustcolor("blue", alpha.f = 0.05)) {
	branch_cols <- grep("^branch_", colnames(df))
	nedge <- nrow(tree$edge)

	branch_color <- c("blue", "red", "green", "orange", "cyan", "brown", "darkgreen", "deeppink")
	branch_color <- rep(branch_color, length.out = nedge)

	bands <- compute_band_limits(0, max(df$time), n_bands)

	add_bands <- function() {
		u <- par("usr")
		for(i in seq_len(nrow(bands))) {
			rect(bands$start[i], u[3], bands$end[i], u[4],
			     col = time_band_col, border = NA)
		}
	}

op <- par(no.readonly = TRUE)
on.exit(par(op))

layout(matrix(1:4, ncol = 1), heights = c(2, 3, 2, 2))

	## Panel 1: tree
	par(mar = c(0, 4, 1, 0.1))

	plot(tree,
	     direction = "rightwards",
	     show.tip.label = FALSE,
	     edge.color = branch_color,
	     edge.width = 3)

	add_bands()
	## Ici les bandes passent légèrement au-dessus de l'arbre.
	## C'est moins élégant, mais beaucoup plus robuste que add = TRUE.

	## Panel 2: branch values
	par(mar = c(0, 4, 0, 0.1))

	plot(df$time, df$branch_1,
	     type = "n",
	     ylim = range(df[, branch_cols]*2, na.rm = TRUE),
	     xlab = "",
	     ylab = "",
	     bty = "n",
	     xaxt = "n",
	     yaxt = "n",
	     cex.lab = 1.5,
	     mgp = c(1, 0.6, 0))

	add_bands()

	mtext("Brownian trait", side = 2, line = 0, cex = 1)
	lines(df$time, df$branch_1,
	      col = adjustcolor(branch_color[1], alpha.f = 0.7),
	      lwd = 1)

	if (nedge >= 2) {
		for (j in 2:nedge) {
			lines(df$time, df[[paste0("branch_", j)]],
			      col = adjustcolor(branch_color[j], alpha.f = 0.7),
			      lwd = 1)
		}
	}

## Panel 3: empirical mean
par(mar = c(0, 4, 0, 0.1))

plot(df$time, df[, "empirical_mean"],
     type = "n",
     ylim = range(df[, "empirical_mean"]*2, na.rm = TRUE),
     xlab = "",
     ylab = "",
     bty = "n",
     cex.lab = 1.75,
     cex.axis = 1.5,
     xaxt = "n",
     yaxt = "n",
     mgp = c(1, 0.6, 0))

add_bands()

mtext("Empirical mean", side = 2, line = 0, cex = 1)
lines(df$time, df[, "empirical_mean"],
      col = "black",
      lwd = 1)

	## Panel 4: empirical variance
	par(mar = c(4, 4, 0, 0.1))

	plot(df$time, df[, "empirical_variance"],
	     type = "n",
	     ylim = range(df[, "empirical_variance"]*2, na.rm = TRUE),
	     xlab = "Time",
	     ylab = "",
	     bty = "n",
	     cex.lab = 1.75,
	     cex.axis = 1.5,
	     xaxt = "s",
	     yaxt = "n",
	     mgp = c(1, 0.6, 0))

	add_bands()
	mtext("Empirical variance", side = 2, line = 0, cex = 1)
	lines(df$time, df[, "empirical_variance"],
	      col = "black",
	      lwd = 1)
}


# Example:
# library(ape)
# library(phytools)
# library(CompQuadForm)
#
# source("bm_empvar_phylo_ext.R")
#
# # Complete tree with extinctions, for example from TreeSim:
# # tr <- TreeSim::sim.bd.age(age = 5, lambda = 1, mu = 0.4,
# #                           numbsim = 1, mrca = TRUE, complete = TRUE)[[1]]
# # sim <- simulate_empvar_bm(tr, n = 200, sigma2 = 1, step = 0.02)
# # plot_empvar_bm(tr, sim)
