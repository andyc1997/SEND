#' @export
eps_estim <- function(X){
  # CHECK number of obs > 1
  n <- dim(X)[1]
  Rk <- rep(0, n-1)

  if (n <= 1) {
    print('n must be greater than 1.')
    return(NULL)
  }

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

    if (is.null(nrow(X1))) X1 <- matrix(X1, nrow = 1)
    for (j in 1:nrow(X1)) {
      r.temp <- min(sqrt(rowSums((sweep(X2, MARGIN = 2, X1[j, ]))^2)))

      if (r.temp < r) {
        r <- r.temp
        js <- j
      }

    }

    Rk[k] <- r
    X2 <- rbind(X2, X1[js, ])
    X1 <- X1[-js,]
  }

  eps <- max(Rk)/2
  return(list(eps = eps, Rk = Rk))
}


#' @export
alarm_estim <- function(newx, X, eps){
  # Alarm raising
  Tn <- min(sqrt(rowSums((sweep(X, MARGIN = 2, newx))^2)))
  Tn <- Tn/eps
  return(Tn)
}


#' @export
sm.boot <- function(X, eps, N.boot, N.sim) {
  # Smoothed Bootstrap method
  alarm.boot <- rep(0, N.boot)

  for (i in 1:N.boot){
    # SAMPLE
    idx.boot <- sample(N.sim, size = N.sim+1, replace = TRUE)

    # SMOOTH
    Z <- uniformly::runif_in_sphere(N.sim+1, d = ncol(X), r = eps)
    sample.boot <- X[idx.boot, ] + Z
    alarm.boot[i] <- alarm_estim(sample.boot[N.sim+1, ], sample.boot[-(N.sim+1), ], eps)
  }

  return(alarm.boot)
}


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
