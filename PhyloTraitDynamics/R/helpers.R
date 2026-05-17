## Internal helpers for PhyloTraitDynamics.
## These functions provide a conservative public interface around the original
## scripts stored in inst/scripts/original/.  The original scientific code is
## sourced in isolated environments and is not edited here.

.bd_original_envs <- new.env(parent = emptyenv())

.bd_original_file <- function(file) {
  path <- system.file("scripts", "original", file,
                      package = "PhyloTraitDynamics", mustWork = FALSE)
  if (!nzchar(path)) {
    ## Fallback useful when the package is used from an unpacked source tree.
    candidates <- c(
      file.path("inst", "scripts", "original", file),
      file.path("..", "inst", "scripts", "original", file)
    )
    ok <- candidates[file.exists(candidates)]
    if (length(ok) > 0L) {
      path <- ok[1L]
    }
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Could not find original script: ", file, call. = FALSE)
  }
  path
}

.bd_source_original <- function(key, file) {
  if (!exists(key, envir = .bd_original_envs, inherits = FALSE)) {
    env <- new.env(parent = globalenv())
    sys.source(.bd_original_file(file), envir = env)
    assign(key, env, envir = .bd_original_envs)
  }
  get(key, envir = .bd_original_envs, inherits = FALSE)
}

.bd_env_empvar <- function() {
  .bd_source_original(
    "empvar",
    "time_dependent_birth_death_brownian_v5_closed_constant.R"
  )
}

.bd_env_meanvar <- function() {
  .bd_source_original(
    "meanvar",
    "time_dependent_birth_death_empirical_mean_variance.R"
  )
}

.bd_env_mrca <- function() {
  .bd_source_original("mrca", "MRCA_age_dynamics.R")
}

.bd_env_fixed <- function() {
  .bd_source_original("fixed", "bm_empvar_phylo_ext.R")
}

.as_rate_function <- function(rate, name = "rate") {
  if (is.numeric(rate) && length(rate) == 1L && is.finite(rate) && rate >= 0) {
    force(rate)
    return(function(t) rep(rate, length(as.numeric(t))))
  }

  if (is.function(rate)) {
    force(rate)
    force(name)
    return(function(t) {
      t <- as.numeric(t)
      y <- try(rate(t), silent = TRUE)

      if (inherits(y, "try-error") || !(length(y) %in% c(1L, length(t)))) {
        y <- try(vapply(t, rate, numeric(1)), silent = TRUE)
      }

      if (inherits(y, "try-error")) {
        stop(name, " must be a numeric constant or a function returning a scalar or one value per input time.", call. = FALSE)
      }

      y <- as.numeric(y)
      if (length(y) == 1L && length(t) > 1L) {
        y <- rep(y, length(t))
      }
      if (length(y) != length(t)) {
        stop(name, " must return either a scalar or one value per input time.", call. = FALSE)
      }
      if (any(!is.finite(y)) || any(y < 0)) {
        stop(name, " must be finite and non-negative on the working time grid.", call. = FALSE)
      }
      y
    })
  }

  stop(name, " must be either a non-negative numeric constant or a function of time.", call. = FALSE)
}

.make_public_grid <- function(time_start = 0, time_end, time_step) {
  if (!is.numeric(time_start) || length(time_start) != 1L ||
      !is.finite(time_start) || time_start < 0) {
    stop("'time_start' must be a finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(time_end) || length(time_end) != 1L ||
      !is.finite(time_end) || time_end < time_start) {
    stop("'time_end' must be a finite number >= 'time_start'.", call. = FALSE)
  }
  if (!is.numeric(time_step) || length(time_step) != 1L ||
      !is.finite(time_step) || time_step <= 0) {
    stop("'time_step' must be a finite strictly positive number.", call. = FALSE)
  }

  grid <- seq(time_start, time_end, by = time_step)
  if (length(grid) == 0L || abs(utils::tail(grid, 1L) - time_end) > 1e-14) {
    grid <- c(grid, time_end)
  }
  unique(grid)
}

.as_source_sim_empvar <- function(x) {
  list(
    grid = x$time,
    sigma2 = x$sigma2,
    V = x$empirical_variance,
    N = x$n_lineages,
    B = x$B,
    tmax = x$time_end,
    dt = x$time_step
  )
}

.as_source_sim_mean <- function(x) {
  list(
    grid = x$time,
    sigma2 = x$sigma2,
    x0 = x$x0,
    M = x$empirical_mean,
    N = x$n_lineages,
    B = x$B,
    tmax = x$time_end,
    dt = x$time_step
  )
}

.match_grid_indices <- function(source_grid, target_grid, tol = 1e-10) {
  idx <- integer(length(target_grid))
  for (i in seq_along(target_grid)) {
    j <- which(abs(source_grid - target_grid[i]) <= tol * max(1, abs(target_grid[i])))
    if (length(j) == 0L) {
      stop(
        "The requested time grid is not available in the simulated grid. ",
        "For simulations, choose 'time_start' and 'time_step' so that requested times lie on seq(0, time_end, by = time_step).",
        call. = FALSE
      )
    }
    idx[i] <- j[1L]
  }
  idx
}

.active_lineage_counts_fixed_tree <- function(tree, times, env) {
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  depths <- env$compute_node_depths(tree)
  edge_start <- depths[tree$edge[, 1]]
  edge_end <- depths[tree$edge[, 2]]
  vapply(times, function(t) sum(edge_start < t & edge_end >= t), integer(1))
}
