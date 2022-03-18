install.packages("readxl")
install.packages("tidyr")
install.packages("readr")
install.packages("sjPlot")
install.packages("sjmisc")
install.packages("ggplot2")
library(readxl)             
library(reshape2)
library(tidyr)
library(dplyr)
library(readr)
library(sjPlot)
library(sjmisc)
library(ggplot2)
#### Pregunta 1
### 1.1 Defina su directorio de trabajo
getwd()

### 1.2 Cargar la data "mRNA_export.csv"
file.choose()
dt <- "C:\\Users\\RICARDO\\Documents\\GitHub\\GBI6_ExamenFinal\\data\\mRNA_expr.csv"
df <- read.csv(dt)

### 1.3 Cree una función de nombre long_df para transformar la data df de la forma wide a long.
long_df <-function(df) {
  library(tidyr)
  library(dplyr)
  df$dataset <- as.factor(df$dataset)
  long_data <- gather(df, gen, expresion_level, GATA3, PTEN, XBP1, ESR1, MUC1, FN1, GAPDH, factor_key=TRUE)
  long_data1 <- select(long_data, -1)
  return (long_data1)
}

### 1.4 Genere la data df_long utilizandola función long_df
df_long <- long_df(df)
test <- df_long %>% group_by(gen)

### 1.5 A partir de la data df_long, genere un Descriptive table (librería sjPlot)
### y guárdelo en el subdirectorio "result" con el nombre mRNA_expr_summary.doc 
### (previamente debe seleccionar las celdas adecuadas con funciones de la librería
### dplyr en conjunto con el operador %>%).
c_m <- test %>% group_by(gen) %>% mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
  select(-row)

descr(c_m) %>% tab_df(title = "Estadistica descriptiva del cancer de MAMA y su nivel de expresion de los diferentes genes",
                              file = "mRNA_expr_summary.doc") 

### 1.6 Interprete los resultados de la tabla resultante de la sección 1.5.


#### Pregunta 2 VISUALIZACION DE DATOS
### 2.1 Cree la funcion tcga_boxplots para visualizar boxplots y jitterplots. 
### El único parámetro de la función es dataset. La función debe crear un boxplot 
### de los niveles de expresión para cada gen.
tcga_boxplots <- function(dataset){
  return(ggplot (dataset, aes(x = gen, y = expresion_level, color = gen)) + labs(x="Gen", y = "Expression level")
         + geom_boxplot() + geom_jitter(size=0.20))
}
tcga_boxplots(df_long)

### 2.2 Usando la función lapply() genera la lista cancertype_boxplots con 
### las gráficas para cada tipo de cancer (cada tipo de dataset).
tipo_de_cancer <- list(
  BRCA <- filter(df_long, dataset == "BRCA"),
  OV <- filter(df_long, dataset == "OV"),
  LUSC <- filter(df_long, dataset == "LUSC"),
  KIPAN <- filter(df_long, dataset == "KIPAN"),
  KIRP <- filter(df_long, dataset == "KIRP"),
  UCEC <- filter(df_long, dataset == "UCEC")
)

cancertype_boxplots = lapply(tipo_de_cancer, function(x) {ggplot(x, aes(dataset, expresion_level, col = dataset)) +
    geom_boxplot() + theme(legend.position='none') + labs(x = "", y = "Expression level")})

### 2.3 Guarde en el subdirectorio result el tercer plot de la lista 
### cancertype_boxplots con el nombre boxplot3.png. La definición de a figura 
### debe ser de 300dpi.
cancertype_boxplots[3]
ggsave("boxplot3.png", width = 6, height = 8, dpi = 300)

### 2.4 Escriba la función reg_gen_expression, con el parámetro gen. En este caso 
### la función visualizará una gráfica de nube de puntos (geom_point) y una 
### regresión por tipo de "dataset". La gráfica será de comparación de gen1 con 
### gen2; por ejemplo en el eje "x" GATA3 y en el eje "y"" FN1.
data_Sett <- list(
  BRCA <- filter(df_long, dataset == "BRCA"),
  OV <- filter(df_long, dataset == "OV"),
  LUSC <- filter(df_long, dataset == "LUSC"),
  KIPAN <- filter(df_long, dataset == "KIPAN"),
  KIRP <- filter(df_long, dataset == "KIRP"),
  UCEC <- filter(df_long, dataset == "UCEC"), 
  GAPDH <- filter(df_long, dataset == "GAPDH")
)

graficasc <- list (  
  BRCA_1 <- BRCA %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row),
  OV_1 <- OV %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row),
  LUSC_1 <- LUSC %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row),
  KIPAN_1 <- KIPAN %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row),
  KIRP_1 <- KIRP %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row),
  UCEC_1 <- UCEC %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row),
  GAPH_1 <- UCEC %>% group_by(gen) %>% mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = gen, values_from = expresion_level) %>%
    select(-row)
)

reg_gen_expression <- function (gen) {
  if (gen == "GATA3"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(GATA3, PTEN)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GATA3, XBP1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GATA3, ESR1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GATA3, MUC1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GATA3, FN1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GATA3, GAPDH)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
  else if (gen == "PTEN"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(PTEN, GATA3)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(PTEN, XBP1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(PTEN, ESR1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(PTEN, MUC1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(PTEN, FN1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(PTEN, GAPDH)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
  else if (gen == "XBP1"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(XBP1, GATA3)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(XBP1, PTEN)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(XBP1, ESR1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(XBP1, MUC1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(XBP1, FN1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(XBP1, GAPDH)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
  else if (gen == "ESR1"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(ESR1, GATA3)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(ESR1, PTEN)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(ESR1, XBP1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(ESR1, MUC1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(ESR1, FN1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(ESR1, GAPDH)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
  else if (gen == "MUC1"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(MUC1, GATA3)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(MUC1, PTEN)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(MUC1, XBP1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(MUC1, ESR1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(MUC1, FN1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(MUC1, GAPDH)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
  else if (gen == "FN1"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(FN1, GATA3)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(FN1, PTEN)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(FN1, XBP1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(FN1, ESR1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(FN1, MUC1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(FN1, GAPDH)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
  else if (gen == "GAPDH"){
    total <- list (
      lapply(graficasc, function(x) {ggplot(x, aes(GAPDH, GATA3)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GAPDH, PTEN)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GAPDH, XBP1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GAPDH, ESR1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GAPDH, MUC1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}), 
      lapply(graficasc, function(x) {ggplot(x, aes(GAPDH, FN1)) + 
          geom_point() + geom_smooth(method='lm', formula= y~x)}))
  }
}
### 2.6. Guarde en el subdirectorio result el onceavo plot de la lista reg_genplots
### con el nombre regresion11.pdf. La definición de a figura debe ser de 300dpi.
graficos_de_cancer <- reg_gen_expression("GATA3")
graficos_cancer[[5]][[2]]
ggsave("regresion11.pdf", width = 6, height = 8, dpi = 300)

### 2.7. Interprete los resultados de las figuras boxplot3.png y regresion11.pdf.