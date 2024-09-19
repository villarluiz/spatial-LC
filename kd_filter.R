# @title Kalman Filter with Drift
#
# @description  Computes mean and variance of the distribution of the state, conditional
# on the covariances of observation and system errors, evolution matrices, drifts
# and the observations up to that point.
#
# @param y  Observed data as an m-by-N matrix where m is the dimension
# @param m0 Initial mean as a p-vector
# @param C0 Initial covariance p-by-p matrix
# @param V  Covariance m-by-p matrix V referring to observation error
# @param W  Covariance p-by-p matrix W referring to system error
# @param Ft Constant observation matrix as an m-by-p matrix
# @param Gt Constant evolution matrix as a p-by-p matrix
# @param a1 Drift of the observational equation
# @param a2 Drift of the system equation
#
# @return A list containing the means `m` and the covariances `C`
#
# @note
# Reference: Petris et al, 2009, p.53
#
kd.filter <- function(y, m0, C0, V, W, Ft, Gt, a1, a2) {
	N <- ncol(y)    # Number of observations
	p <- length(m0) # State dimension

	ret <- vector("list", 2)            # List to be returned
	ret$m <- matrix(nrow = p, ncol = N) # Array of means m_t
	ret$C <- array(dim = c(p, p, N))    # Array of covariances C_t

	# First-step computation

	# at <- a2 + Gt %*% m0
	# Rt <- W + Gt %*% C0 %*% t(Gt)
	# Qt <- Ft %*% Rt %*% t(Ft) + V
	# ybart <- y[ ,1] - a1
	# Bt <- Rt %*% t(Ft)
	# ret$C[ , ,1] <- Rt - Bt %*% solve(Qt) %*% t(Bt)
	# ret$m[ ,1] <- ret$C[ , ,1] %*% (t(Ft) %*% solve(V) %*% ybart + solve(Rt) %*% at)
	at <- a2 + Gt %*% m0
	Rt <- W + Gt %*% C0 %*% t(Gt)
	ft <- Ft %*% at
	Qt <- Ft %*% Rt %*% t(Ft) + V
	et <- y[ ,1] - a1 - ft
	At <- Rt %*% t(Ft) %*% chol2inv(chol(Qt))
	ret$m[ ,1] <- at + At %*% et
	ret$C[ , ,1] <- Rt - At %*% Ft %*% Rt

	# Updating

	# FT.V.inv <- t(Ft) %*% solve(V)
	# ybart <- y - a1

	for (i in 2:N) {
		# at <- a2 + Gt %*% ret$m[ ,i-1]
		# Rt <- W + Gt %*% ret$C[ , ,i-1] %*% t(Gt)
		# Qt <- Ft %*% Rt %*% t(Ft) + V
		# Bt <- Rt %*% t(Ft)
		# ret$C[ , ,i] <- Rt - Bt %*% solve(Qt) %*% t(Bt)
		# ret$m[ ,i] <- ret$C[ , ,i] %*% (FT.V.inv %*% ybart[ ,i] + solve(Rt) %*% at)
		at <- a2 + Gt %*% ret$m[ ,i-1]
		Rt <- W + Gt %*% ret$C[ , ,i-1] %*% t(Gt)
		ft <- Ft %*% at
		Qt <- Ft %*% Rt %*% t(Ft) + V
		et <- y[ ,i] - a1 - ft
		At <- Rt %*% t(Ft) %*% chol2inv(chol(Qt))
		ret$m[ ,i] <- at + At %*% et
		ret$C[ , ,i] <- Rt - At %*% Ft %*% Rt
	}

	ret
}
