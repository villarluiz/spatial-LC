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

data_ob <- read.table("Dados_microrregioes/obt_OPS/obt_2022.txt", sep = ";")

### data evolution 1980-2000-2020-2021 -----
mx <- read.table("Dados_microrregioes/det/logmx_det/logmx_1980.txt")
mx2 <- read.table("Dados_microrregioes/det/logmx_det/logmx_2000.txt")
mx3 <- read.table("Dados_microrregioes/det/logmx_det/logmx_2020.txt")
mx4 <- read.table("Dados_microrregioes/det/logmx_det/logmx_2021.txt")

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

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt1 + plt2 + plt3 + plt4 +
  plot_annotation(title = "Mortality rates for 25-29 age bracket, RJ") +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                                    low = "yellow", high = "red") -> pt
ggsave("mx_25yo_res.png", pt, device = png(width = 800, height = 440))


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

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt1 + plt2 + plt3 + plt4 +
  plot_annotation(title = "Mortality rates for 50-54 age bracket, RJ") +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "yellow", high = "red") -> pt
ggsave("mx_50yo_res.png", pt, device = png(width = 800, height = 440))


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

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt1 + plt2 + plt3 + plt4 +
  plot_annotation(title = "Mortality rates for 80+ age bracket, RJ") +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "yellow", high = "red") -> pt
ggsave("mx_80yo_res.png", pt, device = png(width = 800, height = 440))


### average mortality RJ todo x ano

#### salvar todos os logmor, fazer apply(mean,1 ou 2) e plotar por faixa etaria e por ano. 
avg.year <- c()
aux.age <- c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+")
for(i in 1:42){
  aux <- as.character(i + 1979)
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", aux, ".txt"))[,-1] %>%
    mutate_all(function(x) as.numeric(x)) %>%   apply(1, mean, na.rm = T)
  avg.year <- rbind(avg.year, cbind(mx, aux.age, year = rep(i+1979, 13)))
}
avg.year <- as.data.frame(avg.year)
ggplot(avg.year) + geom_line(aes(x = as.numeric(year), y = as.numeric(mx), color = aux.age)) + labs(x = "year", y = "logmx", color = "age group")
ggsave("graficos/avgmortality_x_year.png", device = png(width = 800, height = 440))
