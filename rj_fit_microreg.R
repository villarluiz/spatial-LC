library(sf)
library(spdep)

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
thin = 5

m0 <- 0; C0 <- 10

### y_obs
y <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:t){
  aux <- i + 1979
  mx <- as.data.frame(read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt")))[,-1] %>%
    stack()
  y[,i] <- mx$values
}
ind <- as.character(mx$ind)

### missing data input (interpolando no ano, nao no espaço)
y_treat <- t(zoo::na.approx(t(y)))
y_treat[92, 42] <- y_treat[92, 41] ## repetindo valor do ultimo ano (2021) pq nao tem como interpolar

# source("sblc_fun.R")
source("sffbs_fun.R")
fit <- sblc(y_treat, ages, t, n, m0, C0, M, W, it = 1000, bn = 500, thin = 1)

### vis. cadeias

plot.ts(fit$kappa.chain[ ,1])
plot.ts(fit$kappa.chain[ ,10])
plot.ts(fit$kappa.chain[ ,20])
plot.ts(fit$kappa.chain[ ,30])

acf(fit$kappa.chain[ ,1])
acf(fit$kappa.chain[ ,10])
acf(fit$kappa.chain[ ,20])
acf(fit$kappa.chain[ ,30])

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
### ok sigma_t ser 3?
