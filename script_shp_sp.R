library(sf)
library(spdep)

rj_micro <- st_read("SP_Microrregioes_2022/SP_Microrregioes_2022.shp")
plot(rj_mun[,c(2,5)])
plot(rj_micro)

shape_aisp_nb_val <- st_make_valid(rj_micro[,c(2,5)])
nb <- poly2nb(shape_aisp_nb_val, queen=TRUE)

# matriz de vizinhaça
W <- nb2mat(nb, style="B", zero.policy = TRUE)

data_ob <- read.table("Dados_microrregioes/obt_OPS/obt_2022.txt", sep = ";")

### data evolution 1980-2000-2020-2021
mx <- read.table("Dados_microrregioes_SP/det/logmx_det/logmx_1980.txt")
mx2 <- read.table("Dados_microrregioes_SP/det/logmx_det/logmx_2000.txt")
mx3 <- read.table("Dados_microrregioes_SP/det/logmx_det/logmx_2020.txt")
mx4 <- read.table("Dados_microrregioes_SP/det/logmx_det/logmx_2021.txt")

##### heatmap 25-29
data.aux <- cbind(rj_micro, t(mx[2,-1]), t(mx2[2,-1]), t(mx3[2,-1]), t(mx4[2,-1]))
data.aux.rang <- cbind(t(mx[2,-1]), t(mx2[2,-1]), t(mx3[2,-1]), t(mx4[2,-1]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2)))) + labs(fill = "mx") +
  ggtitle("1980") -> plt1
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2.1)))) + labs(fill = "mx") +
  ggtitle("2000") -> plt2
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2.2)))) + labs(fill = "mx") +
  ggtitle("2020") -> plt3
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2.3)))) + labs(fill = "mx") +
  ggtitle("2021") -> plt4

lims <- c(min(exp(data.aux.rang), na.rm = T), max(exp(data.aux.rang), na.rm = T))

library(patchwork)
plt1 + plt2 + plt3 + plt4 +
  plot_annotation(title = "Mortality rates for 25-29 age bracket, SP") +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                                    low = "yellow", high = "red") -> pt
ggsave("mx_25yo_res_sp.png", pt, device = png(width = 800, height = 440))


##### heatmap 50-54
data.aux <- cbind(rj_micro, t(mx[7,-1]), t(mx2[7,-1]), t(mx3[7,-1]), t(mx4[7,-1]))
data.aux.rang <- cbind(t(mx[7,-1]), t(mx2[7,-1]), t(mx3[7,-1]), t(mx4[7,-1]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7)))) + labs(fill = "mx") +
  ggtitle("1980") -> plt1
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7.1)))) + labs(fill = "mx") +
  ggtitle("2000") -> plt2
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7.2)))) + labs(fill = "mx") +
  ggtitle("2020") -> plt3
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7.3)))) + labs(fill = "mx") +
  ggtitle("2021") -> plt4

lims <- c(min(exp(data.aux.rang), na.rm = T), max(exp(data.aux.rang), na.rm = T))

library(patchwork)
plt1 + plt2 + plt3 + plt4 +
  plot_annotation(title = "Mortality rates for 50-54 age bracket, SP") +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "yellow", high = "red") -> pt
ggsave("mx_50yo_res_sp.png", pt, device = png(width = 800, height = 440))


##### heatmap 80+
data.aux <- cbind(rj_micro, t(mx[13,-1]), t(mx2[13,-1]), t(mx3[13,-1]), t(mx4[13,-1]))
data.aux.rang <- cbind(t(mx[13,-1]), t(mx2[13,-1]), t(mx3[13,-1]), t(mx4[13,-1]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13)))) + labs(fill = "mx") +
  ggtitle("1980") -> plt1
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13.1)))) + labs(fill = "mx") +
  ggtitle("2000") -> plt2
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13.2)))) + labs(fill = "mx") +
  ggtitle("2020") -> plt3
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13.3)))) + labs(fill = "mx") +
  ggtitle("2021") -> plt4

lims <- c(min(exp(data.aux.rang), na.rm = T), max(exp(data.aux.rang), na.rm = T))

library(patchwork)
plt1 + plt2 + plt3 + plt4 +
  plot_annotation(title = "Mortality rates for 80+ age bracket, SP") +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "yellow", high = "red") -> pt
ggsave("mx_80yo_res_sp.png", pt, device = png(width = 800, height = 440))
