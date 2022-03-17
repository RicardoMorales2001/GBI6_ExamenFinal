install.packages("readxl")
install.packages("tidyr")
library(readxl)             
library(reshape2)
library(tidyr)
library(dplyr)
# Pregunta 1
# 1.1 Defina su directorio de trabajo
getwd()

# 1.2 Cargar la data "mRNA_export.csv"
file.choose()
dt <- "C:\\Users\\RICARDO\\Documents\\GitHub\\GBI6_ExamenFinal\\data\\mRNA_expr.csv"
df <- read.csv(dt)

# 1.3 Cree una función de nombre long_df para transformar la data df de la forma wide a long.
long_data <- gather(df,  GEN, expresion_level,GATA3, PTEN, XBP1, ESR1, MUC1, FN1, GAPDH, factor_key=TRUE)
long_data1 <- select(long_data, -bcr_patient_barcode)
View(long_data1)
