library(tidyverse)

### TESTE
# data <- as.data.frame(t(read.table("Dados_microrregioes/ex_puro/ex_2021.txt", sep = ";")))
# colnames(data) <- c("Idade", as.character(seq(33001, 33018)), "total")
# data <- cbind(data[-c(1,83),-1], "grupo" = c(1, rep(2,4), rep(3,10), rep(4,10), rep(5,10), rep(6,10),
#                            rep(7,10), rep(8,10), rep(9,10), rep(10,6)))
# data <- as.data.frame(sapply(data, as.numeric)) %>%
#   group_by(grupo) %>%
#   summarise_all(sum)
# write.table(data2, "Dados_microrregioes/ex_OPS/ex_2021_OPS.txt")              

### GRUPAMENTO OPS
grupo <- c(1, rep(2,4), rep(3,10), rep(4,10), rep(5,10), rep(6,10),
           rep(7,10), rep(8,10), rep(9,10), rep(10,6))

### FUNCAO
for(i in seq(2021,2011)){
  data <- as.data.frame(t(read.table(paste("Dados_microrregioes/ex_puro/ex_", as.character(i), ".txt", sep = ""), sep = ";")))
  colnames(data) <- c("Idade", as.character(seq(33001, 33018)), "total")
  data <- cbind(data[-c(1,83),-1], "grupo" = grupo)
  data <- as.data.frame(sapply(data, as.numeric)) %>%
    group_by(grupo) %>%
    summarise_all(sum)
  write.table(data, paste("Dados_microrregioes/ex_OPS/ex_", as.character(i), ".txt", sep = ""))
}
