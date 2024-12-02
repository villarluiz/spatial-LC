#' @title Lee-Carter Bayesian Estimation for mortality data
#'
#' @description Performs Bayesian estimation of the Lee-Carter model considering different
#' variances for each age.
#'
#' @usage
#' blc(Y, prior = NULL, init = NULL, M = 5000, bn = 4000, thin = 1)
#'
#' @param Y Log-mortality rates for each age.
#' @param prior A list containing the prior mean \eqn{m_0} and the prior
#' variance \eqn{C_0}.
#' @param init A list containing initial values for \eqn{\alpha}, \eqn{\beta},
#' \eqn{\phi_V}, \eqn{\phi_W} and \eqn{\theta}.
#' @param M The number of iterations. The default value is 5000.
#' @param bn The number of initial iterations from the Gibbs sampler that should be discarded (burn-in). The default value is 4000.
#' @param thin A Positive integer specifying the period for saving samples. The default value is 1.
#'
#' @details
#' Let \eqn{Y_{it}} be the log mortality rate at age \eqn{i} and time \eqn{t}. The Lee-Carter
#' model is specified as follows:
#'
#' \eqn{Y_{it} = \alpha_i + \beta_i \kappa_t + \varepsilon_{it}, i=1,...,p} and \eqn{t=1,...,T},
#'
#' where \eqn{\alpha = (\alpha_1,...,\alpha_p)'} are the interept of the model that represent
#' the log-mortality rate mean in each age; \eqn{\beta = (\beta_1,...,\beta_p)'} are the
#' coefficient regression that represent the speed of relative change in the log-mortality
#' rate in each age. \eqn{\kappa = (\kappa_1,...,\kappa_T)'} are the state variable that
#' represent the global relative change in log-mortality rate. Finally, \eqn{\varepsilon_{it} ~ N(0, \sigma^2_i)}
#' is the random error.
#'
#' For the state variable \eqn{\kappa_t} Lee and Carter (1992) proposed a random walk with
#' drift to govern the dynamics over time:
#'
#' \eqn{\kappa_t = \kappa_{t-1} + \theta + \omega_t},
#'
#' where \eqn{\theta} is the drift parameter and \eqn{\omega_t} is the random error of the
#' random walk.
#'
#' We implement the Bayesian Lee Carter (BLC) model, proposed by Pedroza (2006), to estimate
#' the model. In this approach, we take advantage of the fact that the Bayesian Lee Carter
#' can be specified as dynamic linear model, to estimate the state variables \eqn{\kappa_t}
#' through FFBS algorithm. To estimate the others parameters we use Gibbs sampler to sample
#' from their respective posterior distribution.
#'
#' @return A `BLC` object.
#' \item{alpha}{Posterior sample from alpha parameter.}
#' \item{beta}{Posterior sample from beta parameter.}
#' \item{phiv}{Posterior sample from phiv parameter. phiv is the precision of the random error of the Lee Carter model.}
#' \item{theta}{Posterior sample from theta.}
#' \item{phiw}{Posterior sample from phiw. phiw is the precision of the random error of the random walk.}
#' \item{kappa}{Sampling from the state variables.}
#' \item{Y}{Y Log-mortality rates for each age passed by the user to fit the model.}
#' \item{bn}{The warmup of the algorithm specified by the user to fit the model.}
#' \item{M}{The number of iterations specified by the user to fit the model.}
#' \item{m0}{The prior mean of kappa0.}
#' \item{C0}{The prior covariance matrix of kappa0.}
#'
#' @references Lee, R. D., & Carter, L. R. (1992). “Modeling and forecasting US mortality.” \emph{Journal of the American statistical association}, 87(419), 659-671.
#' @references Pedroza, C. (2006). “A Bayesian forecasting model: predicting US male mortality.” \emph{Biostatistics}, 7(4), 530-550.
#'
#' @examples
#' ## Example of transforming the dataset to fit the function:
#'
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Calculating the mortality rates for the general population:
#' require(dplyr)
#' require(tidyr)
#' require(magrittr)
#'
#' USA %>% mutate(mx = USA$Dx.Total/USA$Ex.Total) -> data
#'
#' data %>%
#' 	filter(Age %in% 18:80 & Year %in% 2000:2015)  %>%
#' 	mutate(logmx = log(mx)) %>%
#' 	dplyr::select(Year,Age,logmx) %>%
#' 	pivot_wider(names_from = Year, values_from = logmx) %>%
#' 	dplyr::select(-Age) %>%
#' 	as.matrix()  %>%
#' 	magrittr::set_rownames(18:80) -> Y
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#' print(fit)
#'
#' ## Viewing the results
#' plot(fit, ages = 18:80)
#' plot(fit, parameter = "beta", ages=18:80)
#' improvement(fit)
#'
#' @import magrittr
#' @import progress
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#'
#' @seealso [fitted.BLC()], [plot.BLC()], [print.BLC()] and [predict.BLC()] for `BLC` methods to native R functions [fitted()],
#'[plot()], [print()] and [predict()].
#'
#'[expectancy.BLC()] and [Heatmap.BLC()] for `BLC` methods to compute and visualise the truncated life expectancy
#'via [expectancy()] and [Heatmap()] functions.
#'
#'[improvement()] to compute the improvement of each age, based on the resulting chains of the beta parameter.
#'
#'
#' @export
blc_missing <- function(Y, prior = NULL, init = NULL, M = 5000, bn = 4000, thin = 1) {
  # -------- Type validation --------
  
  if (mode(Y) != "numeric")
    stop("Expected `Y` to be numeric")
  
  if(!is.null(prior)){
    if (mode(prior) != "list") {stop("Expected `prior` and `init` to be lists")}
  }else{
    prior <- list(m0 = 0, C0 = 100)
  }
  
  if(!is.null(init)){
    if (mode(init) != "list") {stop("Expected `prior` and `init` to be lists")}
  }else{
    init <- list(alpha = runif(nrow(Y)), beta = runif(nrow(Y)), phiv = rep(1, nrow(Y)), phiw = 1, theta = runif(1))
  }
  
  if (mode(prior) != "list" || mode(init) != "list")
    stop("Expected `prior` and `init` to be lists")
  
  
  prior.types <- unique(sapply(prior, mode))
  if (length(prior.types) != 1 || prior.types[1] != "numeric")
    stop("Expected `prior` to only contain numerics")
  
  init.types <- unique(sapply(init, mode))
  if (length(init.types) != 1 || init.types[1] != "numeric")
    stop("Expected `init` to only contain numerics")
  
  if (mode(M) != "numeric" || round(M) != M || M <= 0)
    stop("Expected `M` to be a positive integer")
  
  # if (mode(std.type) != "character" || !(std.type %in% c("incl", "beta")))
  # 	stop("Expected `std` to be one of 'incl' or 'beta'")
  
  if (!(mode(bn) %in% c("NULL", "numeric")))
    stop("Expected `bn` to be either nil or numeric")
  
  # -------- List validation --------
  
  prior.names <- c("C0", "m0")
  init.names <- c("alpha", "beta", "phiv", "phiw", "theta")
  
  if (any(sort(names(prior)) != prior.names))
    stop("Invalid names in argument `prior`")
  
  if (any(sort(names(init)) != init.names))
    stop("Invalid names in argument `init`")
  
  # -------- Dimension validation --------
  
  N <- ncol(Y)
  m <- nrow(Y)
  
  if (length(prior$m0) != 1)
    stop("Invalid dimensions for `prior$m0`")
  
  if (length(prior$C0) != 1)
    stop("Invalid dimensions for `prior$C0`")
  
  if (length(init$alpha) != m)
    stop("Invalid dimensions for `init$alpha`")
  
  if (length(init$beta) != m)
    stop("Invalid dimensions for `init$beta`")
  
  if (length(init$phiv) != m)
    stop("Invalid dimensions for `init$phiv`")
  
  if (length(init$phiw) != 1)
    stop("Invalid dimensions for `init$phiw`")
  
  if (length(init$theta) != 1)
    stop("Invalid dimensions for `init$theta`")
  
  # -------- Initialization --------
  
  std.type = "incl"
  
  ## Progress Bar
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = M, clear = FALSE, width = 60)
  
  # Allocate storage
  ind_missing = which(is.na(Y), arr.ind = T)
  
  chain <- list()
  chain$alpha <- matrix(nrow = m, ncol = M)
  chain$beta <- matrix(nrow = m, ncol = M)
  chain$phiv <- matrix(nrow = m, ncol = M)
  chain$theta <- numeric(M)
  chain$phiw <- numeric(M)
  chain$kappa <- matrix(nrow = N, ncol = M)
  chain$input <- matrix(nrow = nrow(ind_missing), ncol = M)
  
  chain$mt <- matrix(nrow = N, ncol = M)
  chain$Ct <- matrix(nrow = N, ncol = M)
  
  pb$tick() ## Progress Bar
  
  # Initialize chains
  chain$alpha[ ,1] <- init$alpha
  chain$beta[ ,1] <- init$beta
  chain$phiv[ ,1] <- init$phiv
  chain$theta[1] <- init$theta
  chain$phiw[1] <- init$phiw
  
  # Missing treatment
  y_obs = Y
  #
  
  chain$input[,1] <- rowMeans(Y, na.rm = T)[ind_missing[,1]]
  y_obs[ind_missing] <- chain$input[,1]
  
  # -------- Gibbs --------
  
  for (i in 2:M) {
    
    pb$tick() ## Progress Bar
    
    kdf <- kd.filter(y_obs, prior$m0, prior$C0, diag(1/chain$phiv[ ,i-1]),
                     1/chain$phiw[i-1], chain$beta[ ,i-1], 1,
                     chain$alpha[ ,i-1], chain$theta[i-1])
    
    chain$mt[, i] <- kdf$m
    chain$Ct[, i] <- kdf$C
    
    kds <- kd.smoother(y_obs, kdf, 1/chain$phiw[i-1], 1, chain$theta[i-1])
    chain$kappa[ ,i] <- rnorm(N, kds$s[1, ], sqrt(kds$S[1,1, ]))
    
    # Standardize the values
    inclination <- (chain$kappa[1,i] - chain$kappa[N,i]) / (N - 1)
    level <- mean(chain$kappa[ ,i])
    chain$kappa[ ,i] <- (chain$kappa[ ,i] - level) / inclination
    chain$alpha[ ,i-1] <- chain$alpha[ ,i-1] + chain$beta[ ,i-1] * level
    chain$beta[ ,i-1] <- chain$beta[ ,i-1] * inclination
    
    for (j in 1:m) {
      # Generate phiv
      B <- sum((y_obs[j, ] - chain$alpha[j,i-1] - chain$beta[j,i-1] * chain$kappa[ ,i])^2)
      chain$phiv[j,i] <- rgamma(1, N/2, B/2)
      
      # Generate alpha and beta
      X <- cbind(1, chain$kappa[ ,i])
      aux.reg <- chol2inv(chol(t(X) %*% X))
      mean.reg <- aux.reg %*% t(X) %*% y_obs[j, ]
      var.reg <- (1 / chain$phiv[j,i]) * aux.reg
      tmp <- mvrnorm(1, mean.reg, var.reg)
      chain$alpha[j,i] <- tmp[1]
      chain$beta[j,i] <- tmp[2]
    }
    
    # Generate phiw
    B <- sum((chain$kappa[2:N,i] - chain$theta[i-1] - chain$kappa[1:(N-1),i])^2)
    chain$phiw[i] <- rgamma(1, (N-1)/2, B/2)
    
    # Generate theta
    B <- 1 / (chain$phiw[i] * (N-1))
    A <- chain$phiw[i] * (chain$kappa[N,i] - chain$kappa[1,i])
    chain$theta[i] <- rnorm(1, B * A, sqrt(B))
    
    # Missing
    for(j in 1:nrow(ind_missing)){
      loc <- c(ind_missing[j,1],
               ind_missing[j,2])
      
      aux = chain$alpha[loc[1], i] + (chain$beta[loc[1], i] * chain$kappa[loc[2], i])
      chain$input[j, i] <- rnorm(1, mean = aux, sd = sqrt(1/chain$phiv[loc[1], i]))
    }
    y_obs[ind_missing] <- chain$input[,i]

  }
  
  # -------- Return --------
  
  chain$kappa[,1] <- chain$kappa[,2]
  
  chain = lapply(chain, function(x){if(length(dim(x)) > 1){x[,seq(bn+1, M, by = thin)]}else{x[seq(bn+1, M, by = thin)]}})
  
  class(chain) <- "BLC"
  chain$Y <- Y
  chain$bn <- 0
  chain$M <- length(chain$theta) ## Final sample size
  chain$m0 <- prior$m0
  chain$C0 <- prior$C0
  
  chain
}