install.packages("readxl")
install.packages("tidyr")
install.packages("readr")
install.packages("sjPlot")
install.packages("sjmisc")
library(readxl)             
library(reshape2)
library(tidyr)
library(dplyr)
library(readr)
library(sjPlot)
library(sjmisc)
# Pregunta 1
# 1.1 Defina su directorio de trabajo
getwd()

# 1.2 Cargar la data "mRNA_export.csv"
file.choose()
dt <- "C:\\Users\\RICARDO\\Documents\\GitHub\\GBI6_ExamenFinal\\data\\mRNA_expr.csv"
df <- read.csv(dt)

# 1.3 Cree una función de nombre long_df para transformar la data df de la forma wide a long.
long_df <-function(df) {
  library(tidyr)
  library(dplyr)
  df$dataset <- as.factor(df$dataset)
  long_data <- gather(df, gen, expresion_level, GATA3, PTEN, XBP1, ESR1, MUC1, FN1, GAPDH, factor_key=TRUE)
  long_data1 <- select(long_data, -1)
  return (long_data1)
}
View(long_data1)

# 1.4 Genere la data df_long utilizandola función long_df
df_long <- long_df(df)
test <- df_long %>% group_by(gen)

