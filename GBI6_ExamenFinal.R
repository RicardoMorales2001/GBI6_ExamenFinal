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
View(long_data1)

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

descr(c_m) %>% tab_df(title = "Estadistica del cancer de MAMA y el nivel de expresion de los diferentes genes",
                              file = "mRNA_expr_summary.doc") 

### 1.6 Interprete los resultados de la tabla resultante de la sección 1.5.


#### Pregunta 2 VISUALIZACION DE DATOS
### 2.1 Cree la funcion tcga_boxplots para visualizar boxplots y jitterplots. 
### El único parámetro de la función es dataset. La función debe crear un boxplot 
### de los niveles de expresión para cada gen.
tcga_boxplots <- function(dataset){
  library(ggplot2)
  return(ggplot (dataset, aes(x = gen, y = expresion_level, col = gen)) + labs(x="Gen", y = "Expression level")
         + geom_boxplot() + geom_jitter(size=0.20))
}
tcga_boxplots(df_long)
View(tcga_boxplots)