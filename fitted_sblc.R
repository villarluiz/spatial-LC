
fitted_sblc <- function(object, prob = 0.95, ...) {
  obj = object
	N <- dim(obj$alpha.chain)[1]
	t <- obj$info$t
	ages <- obj$info$ages
	n <- obj$info$n
	fits <- array(dim = c(ages*n, t, N))
	

	for (i in 1:N) {
	  A <- matrix(kronecker(rep(1, n), obj$alpha.chain[i,]))
	  B <- matrix(kronecker(rep(1, n), obj$beta.chain[i,]))
	  C <- matrix(kronecker(obj$theta.chain[i,], obj$gamma.chain[i,]))
	  
		fits[ , , i] <- A[, rep(1, t)] + (B %*% obj$kappa.chain[i,]) + C[, rep(1, t)]
	}

	alpha <- 1 - prob
	mean <- apply(fits, c(1,2), mean)
	upper <- apply(fits, c(1,2), quantile, 1 - alpha/2)
	lower <- apply(fits, c(1,2), quantile, alpha/2)

	list(mean = mean, upper = upper, lower = lower)
}



#' fitted.PredBLC <- function(object, prob = 0.95, ...) {
#'   obj = object
#'   fits <- obj$y
#' 
#' 
#'   alpha <- 1 - prob
#'   mean <- apply(fits, c(3,2), mean)
#'   upper <- apply(fits, c(3,2), quantile, 1 - alpha/2)
#'   lower <- apply(fits, c(3,2), quantile, alpha/2)
#' 
#'   colnames(mean) <- colnames(obj$y)
#'   colnames(upper) <- colnames(obj$y)
#'   colnames(lower) <- colnames(obj$y)
#' 
#'   row.names(mean) <- row.names(obj$y)
#'   row.names(upper) <- row.names(obj$y)
#'   row.names(lower) <- row.names(obj$y)
#' 
#'   list(mean = 1 - exp(-exp(mean)), upper = 1 - exp(-exp(upper)), lower = 1 - exp(-exp(lower)))
#' }
