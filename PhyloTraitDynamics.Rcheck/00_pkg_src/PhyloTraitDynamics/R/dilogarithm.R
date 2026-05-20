## ================================================================
## Shared internal utilities for real dilogarithm computations.
##
## These functions are intentionally non-exported. They collect the
## numerical Li_2 routines shared by the empirical-mean and
## empirical-variance modules.
## ================================================================

.dilogarithm_li2_series_0_05 <- function(x, tol = 1e-12, max_iter = 10000L){
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

.dilogarithm_li2_01_scalar <- function(z, tol = 1e-12, max_iter = 10000L){
  if (!is.finite(z)){
    stop("Non-finite argument passed to .dilogarithm_li2_01().")
  }
  if (z < -1e-12 || z > 1 + 1e-12){
    stop(".dilogarithm_li2_01() is implemented only for 0 <= z <= 1.")
  }

  z <- min(1, max(0, z))

  if (z == 0){
    return(0)
  }
  if (z == 1){
    return(pi^2 / 6)
  }
  if (z <= 0.5){
    return(.dilogarithm_li2_series_0_05(z, tol = tol, max_iter = max_iter))
  }

  ## Li_2(z) = pi^2/6 - log(z) log(1-z) - Li_2(1-z)
  pi^2 / 6 - log(z) * log1p(-z) -
    .dilogarithm_li2_series_0_05(1 - z, tol = tol, max_iter = max_iter)
}

.dilogarithm_li2_01 <- function(z, tol = 1e-12, max_iter = 10000L){
  vapply(
    as.numeric(z),
    .dilogarithm_li2_01_scalar,
    numeric(1),
    tol = tol,
    max_iter = as.integer(max_iter)
  )
}

.dilogarithm_li2_real_scalar <- function(z, tol = 1e-12, max_iter = 10000L){
  ## Real dilogarithm Li_2(z) on the real domain z <= 1.
  ## This is enough for the constant birth-death closed forms below.
  if (!is.finite(z)){
    stop("Non-finite argument passed to .dilogarithm_li2_real().")
  }
  if (z > 1 + 1e-12){
    stop(".dilogarithm_li2_real() is implemented only for real arguments <= 1.")
  }
  if (z > 1){
    z <- 1
  }
  if (z == 1){
    return(pi^2 / 6)
  }
  if (z == 0){
    return(0)
  }

  if (z < 0){
    ## Li_2(z) = -Li_2(z/(z-1)) - 1/2 log(1-z)^2.
    w <- z / (z - 1)
    return(-.dilogarithm_li2_real_scalar(w, tol = tol, max_iter = max_iter) -
             0.5 * log1p(-z)^2)
  }

  if (z <= 0.5){
    return(.dilogarithm_li2_series_0_05(z, tol = tol, max_iter = max_iter))
  }

  ## Li_2(z) = pi^2/6 - log(z) log(1-z) - Li_2(1-z).
  pi^2 / 6 - log(z) * log1p(-z) -
    .dilogarithm_li2_real_scalar(1 - z, tol = tol, max_iter = max_iter)
}

.dilogarithm_li2_real <- function(z, tol = 1e-12, max_iter = 10000L){
  vapply(
    as.numeric(z),
    .dilogarithm_li2_real_scalar,
    numeric(1),
    tol = tol,
    max_iter = as.integer(max_iter)
  )
}
