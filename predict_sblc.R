predict.SBLC <- function(object, h, ...) {
  obj = object
	t <- ncol(obj$kappa.chain)  ### N
	ages <- ncol(obj$beta.chain)  #### q
	n <- ncol(obj$theta.chain)
	ss <- (obj$info$it - obj$info$bn)/obj$info$thin

	sim <- array(dim = c(ss, h, n*ages))
	
	pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = ss+1, clear = FALSE, width = 60)
	pb$tick()

	for (l in (1:ss)+1) {
	  
	  pb$tick()
	  
		#fit <- ff_sp(obj$info$y, obj$info$ages, obj$info$t, obj$info$n, obj$info$m0,
		#              obj$info$C0, obj$sigma_e.chain[l,], obj$sigma_w.chain[l],
		#              obj$nu.chain[l], obj$alpha.chain[l,], obj$beta.chain[l,],
		#              obj$gamma.chain[l,], obj$theta.chain[l,])
		
		A <- as.matrix(kronecker(rep(1, n), obj$alpha.chain[l,]))
		B <- as.matrix(kronecker(rep(1, n), obj$beta.chain[l,]))
		C <- as.matrix(kronecker(obj$theta.chain[l,], obj$gamma.chain[l,]))

		at <- obj$other$mt[l, t] + obj$nu.chain[l]
		Rt <- obj$other$Ct[l, t] + obj$sigma_w.chain[l]
		
		ft <- A + (B * at) + C
		Qt <-  (Rt * tcrossprod(B)) + diag(obj$sigma_e.chain[l,], nrow = n*ages)
		
		sim[l-1, 1, ] <- MASS::mvrnorm(1, ft, Qt)

		if (h > 1) for (k in 2:h) {
			at <- at + obj$nu.chain[l]
			Rt <- Rt + obj$sigma_w.chain[l]
			
			ft <- A + (B * at) + C
			Qt <-  (Rt * tcrossprod(B)) + diag(obj$sigma_e.chain[l,], nrow = n*ages)

			sim[l-1, k, ] <- MASS::mvrnorm(1, ft, Qt)
		}
	}

	sim <- list(y = sim, h = h)
	#class(sim) <- "PredBLC"
	sim
}
