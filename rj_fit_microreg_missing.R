library(sf)
library(spdep)
library(tidyverse)

#rj_mun <- st_read("RJ_Municipios_2022/RJ_Municipios_2022.shp")
rj_micro <- st_read("RJ_Microrregioes_2022/RJ_Microrregioes_2022.shp")
#plot(rj_mun[,c(2,5)])
plot(rj_micro[,c(2,5)])

shape_aisp_nb_val <- st_make_valid(rj_micro[,c(2,5)])
nb <- poly2nb(shape_aisp_nb_val, queen=TRUE)

# matriz de vizinhaça
W <- nb2mat(nb, style="B", zero.policy = TRUE)
M <- diag(rowSums(W))

ages <- 13  ## idades    20-24, 25-29, ..., 80+
t <- 42  ## tempo
n <- 18   ## regioes

it = 20000 #iterations
bn = 10000
thin = 4

m0 <- 0; C0 <- 100

### y_obs
y <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:t){
  aux <- i + 1979
  mx <- as.data.frame(read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt")))[,-1] %>%
    stack()
  y[,i] <- mx$values
}
ind <- as.character(mx$ind)

###input missing control
y.control <- y
y[3,3] <- NA; y[12, 12] <- NA

source("sblc_fun.R")
source("sffbs_fun.R")

fit <- sblc_missing(y, ages, t, n, m0, C0, M, W, it = 2000, bn = 1000, thin = 1)
saveRDS(fit, "rj_micro_fit_missing.RDS")

fit <- readRDS("rj_micro_fit_missing.RDS")


plot.ts(fit$alpha.chain[ ,1])
plot.ts(fit$alpha.chain[ ,10])
plot.ts(fit$alpha.chain[ ,13])

acf(fit$alpha.chain[ ,1])
acf(fit$alpha.chain[ ,10])
acf(fit$alpha.chain[ ,13])

plot.ts(fit$beta.chain[ ,1])
plot.ts(fit$beta.chain[ ,10])
plot.ts(fit$beta.chain[ ,13])

acf(fit$beta.chain[ ,1])
acf(fit$beta.chain[ ,10])
acf(fit$beta.chain[ ,13])

plot.ts(fit$gamma.chain[ ,1])
plot.ts(fit$gamma.chain[ ,10])
plot.ts(fit$gamma.chain[ ,13])

acf(fit$gamma.chain[ ,1])
acf(fit$gamma.chain[ ,10])
acf(fit$gamma.chain[ ,13])

plot.ts(fit$theta.chain[ ,1])
plot.ts(fit$theta.chain[ ,10])
plot.ts(fit$theta.chain[ ,13])

acf(fit$theta.chain[ ,1])
acf(fit$theta.chain[ ,10])
acf(fit$theta.chain[ ,13])
graphics.off()


plot.ts(fit$nu.chain)
acf(fit$nu.chain)
plot.ts(fit$sigma_w.chain)
acf(fit$sigma_w.chain)
plot.ts(fit$sigma_e.chain)
acf(fit$sigma_e.chain)
plot.ts(fit$sigma_t.chain)
acf(fit$sigma_t.chain)

plot.ts(apply(fit$theta.chain, 2, median))
### ok sigma_t ser 5?

plot(apply(fit$kappa.chain, 2, median), type = "l")
lines(apply(fit$kappa.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$kappa.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$alpha.chain, 2, median), type = "l")
lines(apply(fit$alpha.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$alpha.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$beta.chain, 2, median), type = "l")
lines(apply(fit$beta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$beta.chain, 2, quantile, 0.975), col = 2, lty = 2)
abline(h = 0)

plot(apply(fit$gamma.chain, 2, median), type = "l")
lines(apply(fit$gamma.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$gamma.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$theta.chain, 2, median), type = "l")
lines(apply(fit$theta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$theta.chain, 2, quantile, 0.975), col = 2, lty = 2)


which(is.na(y), arr.ind = T)

plot.ts(fit$input.chain[,1])
abline(h = y.control[3, 3], col = "red")

plot.ts(fit$input.chain[,3])
abline(h = y.control[12, 12], col = "red")

apply(fit$input.chain, 2, quantile, c(0.025, 0.5, 0.975))[,c(1,3)] -> input_sblc



#### MARGINAL
source("kd_filter.R")
source("kd_smoother.R")
source("blc_missing_mar.R")
Y = y[1:13,]

library(MASS)
library(BayesMortalityPlus)
fit_mar <- blc_missing(Y, M = 2000, bn = 1000)
plot(fit_mar)

plot.ts(fit$input.chain[,1])
lines(fit_mar$input[1,], col = "blue")
abline(h = y.control[3, 3], col = "red")

c("SD_SBLC" = sd(fit$input.chain[,1]), "SD_BLC" = sd(fit_mar$input[1,]))

plot.ts(fit$input.chain[,3])
lines(fit_mar$input[2,], col = "blue")
abline(h = y.control[12, 12], col = "red")

c("SD_SBLC" = sd(fit$input.chain[,3]), "SD_BLC" = sd(fit_mar$input[2,]))

### os erros do modelo blc sao maiores mas o phiv é diferente, testar marginal com phiV
### unico ou sblc com phiv evoluindo na idade

##comparar os qx_fitted

##errobar
library(ggplot2)
# df2 <- 
# ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) +
#   geom_line() +
#   geom_point()+
#   geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
#                 position=position_dodge(0.05))