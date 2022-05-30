#' Estimate eps from the data
#'
#' @param X data matrix contains n observations and p columns.
#' @return The smallest radius eps such that the dataset is connected in the topological space.
#' @export
eps_estim <- function(X){
  # CHECK number of obs > 1
  if (is.null(dim(X)))
    stop('X should contain > 1 observation.')

  n <- dim(X)[1]
  Rk <- rep(0, n-1)

  # SHUFFLE data
  # X1 contains remaining data
  # X2 contains previously selected data
  X1 <- X[sample(n), ]
  X2 <- t(as.matrix(X1[1, ]))
  X1 <- X1[-1, ]

  # COMPUTE eps
  for (k in 1:(n-1)) {
    r <- Inf
    js <- 0

    # vector2matrix
    if (is.null(nrow(X1))) X1 <- matrix(X1, nrow = 1)

    # find min distance
    for (j in 1:nrow(X1)) {
      r.temp <- min(sqrt(rowSums((sweep(X2, MARGIN = 2, X1[j, ]))^2)))

      if (r.temp < r) {
        r <- r.temp
        js <- j
      }
    }

    # update param
    Rk[k] <- r
    X2 <- rbind(X2, X1[js, ])
    X1 <- X1[-js,]
  }

  # final estimate
  eps <- max(Rk)/2
  return(list(eps = eps, Rk = Rk))
}



#' Estimate Tn
#'
#' @param newx a vector for the new observation.
#' @param X data matrix contains n observations and p columns.
#' @param eps a positive scalar which defines the neighborhood of each data point.
#' @return The ratio of Tn to eps
#' @export
alarm_estim <- function(newx, X, eps){
  # CHECK
  if (eps <= 0) stop('eps should be positive.')
  if (is.null(dim(X)))
    stop('X should contain > 1 observation.')

  # Estimator T_n
  Tn <- min(sqrt(rowSums((sweep(X, MARGIN = 2, newx))^2)))
  Tn <- Tn/eps
  return(Tn)
}



#' Perform smoothed boostrap method to estimate the threshold
#'
#' @param X data matrix contains n observations and p columns.
#' @param eps a positive scalar which defines the neighborhood of each data point.
#' @param N.boot the number of bootstrap samples.
#' @export
sm.boot <- function(X, eps, N.boot = 10000) {
  # CHECK
  if (eps <= 0) stop('eps should be positive.')
  if (is.null(dim(X)))
    stop('X should contain > 1 observation.')

  # Smoothed Bootstrap method
  n <- nrow(X)
  alarm.boot <- rep(0, N.boot)

  for (i in 1:N.boot){
    # SAMPLE
    idx.boot <- sample(n, size = n+1, replace = TRUE)

    # SMOOTH
    Z <- uniformly::runif_in_sphere(n+1, d = ncol(X), r = eps)
    sample.boot <- X[idx.boot, ] + Z
    alarm.boot[i] <- alarm_estim(sample.boot[n+1, ], sample.boot[-(n+1), ], eps)
  }

  return(alarm.boot)
}



#' Perform cross-validation smoothing to estimate the threshold
#'
#' @param X data matrix contains n observations and p columns.
#' @param lo the smallest value for eps in the grid.
#' @param hi the largest value for eps in the grid.
#' @param len the number of points considered in the grid.
#' @param alpha a vector of quantiles for eps.
#' @param eps.data a positive scalar which defines the neighborhood of each data point.
#' @return the list of chosen thresholds and the empirical probabilities.
#' @export
cv.smooth <- function(X, lo, hi, len, alpha, eps.data){
  # Cross-validation Smoothing
  # define eps grid
  N.sim <- nrow(X)
  eps.control <- list(lo = lo, hi = hi, len = len)
  eps.grid <- with(eps.control, seq(lo, hi, length.out = len))
  P.hat <- rep(0, eps.control$len)

  # evaluate leave-one-out Tn
  for (k in 1:eps.control$len){
    cnt <- 0

    for (i in 1:N.sim) {
      Tn <- alarm_estim(X[i, ], X[-i, ], eps.grid[k])
      if (Tn > 1) cnt <- cnt + 1
    }

    P.hat[k] <- 1-cnt/N.sim
  }

  # create lookup table
  lookup <- cbind(eps.grid, P.hat)
  eps.min <- sapply(alpha, function(x) lookup[which.min(abs(P.hat - 1 + x)), ])
  return(list(alarm.cv = eps.min[1, ]/eps.data, lookup = lookup))
}
