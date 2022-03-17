install.packages("readxl")
library(readxl)             
library(reshape2)
# Pregunta 1
# 1.1 Defina su directorio de trabajo
getwd()

# 1.2 Cargar la data "mRNA_export.csv"
file.choose()
dt <- "C:\\Users\\RICARDO\\Documents\\GitHub\\GBI6_ExamenFinal\\data\\mRNA_expr.csv"
rd <- read.csv(dt)

# 1.3 Cree una función de nombre long_df para transformar la data df de la forma wide a long.
longrd <- melt(rd, id.vars= "bcr_patient_barcode")  
view(longrd)
