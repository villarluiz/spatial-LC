#' @name predict.BLC
#' @rdname predict.BLC
#'
#' @title BLC: Forecasting
#'
#' @description Calculates the means and variances of the forecast distributions based on
#' the resulting chains from an estimation method.
#'
#'
#' @param object A `BLC` object that is result of a call to blc() function.
#' @param h The prediction horizon.
#' @param ... Other arguments.
#'
#' @return A `PredBLC` object that contains a list with predicted values calculated
#' from `BLC` object chains structured in an array.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' ## Prediction for 2 years ahead
#' pred = predict(fit, h = 2)
#' print(pred)
#'
#' @importFrom MASS mvrnorm
#'
#' @seealso [fitted.BLC()], [print.BLC()], and [plot.PredBLC()] for `PredBLC` methods to native R functions [fitted()],
#'[print()], and [plot()].
#'
#'[expectancy.BLC()] and [Heatmap.BLC] to compute and plot the life expectancy of the prediction(s).
#'
#' @export
predict_missing_BLC <- function(object, h, ...) {
  obj = object
	N <- ncol(obj$Y)
	q <- nrow(obj$beta)

	sim <- array(dim = c(obj$M - obj$bn, h, q))

	for (l in 1:(obj$M - obj$bn)) {
		est.alpha <- obj$alpha[ ,l + obj$bn]
		est.beta <- obj$beta[ ,l + obj$bn]
		est.V <- diag(1/obj$phiv[ ,l + obj$bn])
		est.theta <- obj$theta[l + obj$bn]
		est.W <- 1/obj$phiw[l + obj$bn]

		Gt <- 1

		#filt <- kd.filter(obj$Y, obj$m0, obj$C0, est.V, est.W, est.beta,
		#				  Gt, est.alpha, est.theta)

		a <- est.theta + Gt %*% object$mt[N, l]
		R <- Gt %*% object$Ct[N, l] %*% t(Gt) + est.W
		f <- est.alpha + est.beta %*% as.matrix(a)
		Q <- est.beta %*% as.matrix(R) %*% t(est.beta) + est.V
		sim[l,1, ] <- mvrnorm(1, f, Q)

		if (h > 1) for (k in 2:h) {
			a <- est.theta + Gt %*% a
			R <- Gt %*% R %*% t(Gt) + est.W
			f <- est.alpha + est.beta %*% as.matrix(a)
			Q <- est.beta %*% as.matrix(R) %*% t(est.beta) + est.V

			sim[l,k, ] <- mvrnorm(1, f, Q)
		}
	}

	sim <- list(y = sim, h = h)
	class(sim) <- "PredBLC"
	sim
}
