# @title Kalman Smoother with Drift
#
# @description Computes mean and variance of the distribution of the state, conditional
# on the covariances of observation and system errors, equation matrices
# and all of the observations.
#
# @param y  Observed data as an m-by-N matrix where m is the dimension
# @param f  Kalman filtering results as returned by `k.filter` on the data
# @param W  Covariance p-by-p matrix W referring to system error
# @param Gt Evolution matrix as a p-by-p matrix
# @param a2 Drift of the system equation
#
# @return A list containing the means `s` and the covariances `S`
#
# @note
# Reference: Petris et al, 2009, p.61
#
kd.smoother <- function(y, f, W, Gt, a2) {
	N <- ncol(y)   # Number of observations
	p <- nrow(f$m) # State dimension

	ret <- vector("list", 2)            # List to be returned
	ret$s <- matrix(nrow = p, ncol = N) # Array of means s_t
	ret$S <- array(dim = c(p, p, N))    # Array of covariances S_t

	ret$s[ ,N] <- f$m[ ,N]
	ret$S[ , ,N] <- f$C[ , ,N]

	W.inv <- chol2inv(chol(W))
	GT.W.inv.G <- t(Gt) %*% W.inv %*% Gt
	# G.W.inv <- Gt %*% W.inv
	GT.W.inv <- t(Gt) %*% W.inv
	for (i in (N-1):1) {
		C.inv <- chol2inv(chol(f$C[ , ,i]))
		Bt <- chol2inv(chol(GT.W.inv.G + C.inv))
		ret$s[ ,i] <- Bt %*% (GT.W.inv %*% (ret$s[ ,i+1] - a2) + C.inv %*% f$m[ ,i])
		A <- Bt %*% GT.W.inv
		ret$S[ , ,i] <- Bt + A %*% ret$S[ , ,i+1] %*% t(A)
	}

	ret
}
