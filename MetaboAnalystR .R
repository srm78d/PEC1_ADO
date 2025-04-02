# MetaboAnalyst.R

# Carga de Librerías

library("MetaboAnalystR")
library(vegan)
library(factoextra)
library(pheatmap)
library(scales)

# Establecimiento del objeto 'mSet', estructura central para el manejo de datos y resultados del análisis.
mSet<-InitDataObjects("pktable", "stat", FALSE)

# Importación de datos sin procesar desde el archivo especificado.
mSet<-Read.TextData(mSet, "microbiota_MetaboAnalyst_.txt", "colu", "disc")

# verificar la integridad de los datos importados. Asegurar que los datos estén en un formato adecuado para análisis posteriores y para detectar posibles errores antes de cualquier análisis.
getmSet<-SanityCheckData(mSet)

# Corrección de valores mínimos en el conjunto de datos, paso crucial para prevenir errores en transformaciones matemáticas subsecuentes
mSet<-ReplaceMin(mSet);

# # Selección de variables (metabolitos) según criterios definidos, con el objetivo de minimizar el ruido y destacar elementos significativos en el análisis.Aunque este filtrado no es estrictamente necesario.
mSet<-FilterVariable(mSet, "F", 25, "none", -1, "mean", 0)

# # Adecuación de los datos para la normalización, proceso que optimiza el formato de los mismos
mSet<-PreparePrenormData(mSet)

# # Ajuste de los datos para eliminar variaciones no biológicas entre muestras, mejorando la precisión de las comparaciones.
mSet<-Normalization(mSet, "SumNorm", "NULL", "ParetoNorm", ratio=FALSE, ratioNum=20)

# Genera graficas para evaluar la normalización aplicada.
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)

# Genera graficas para evaluar la normalización de las muestras.
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# realiza un Análisis de Componentes Principales (PCA), una técnica para exploración de la estructura general de los datos y detección de patrones.
mSet<-PCA.Anal(mSet)

# Genera un grafico de pares de componentes principales para visualizar las relaciones entre los componentes principales.
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)

# Genera un grafico de scree plot para determinar cuantos componentes principales son relevantes.
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)

# Genera un grafico 2d de los componentes principales.
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0, "na")

# genera un mapa de calor para visualizar patrones de expresión de variables entre muestras.
mSet<-PlotStaticHeatMap(mSet, "heatmap_1_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, "overview", T, T, NULL, T, F, T, T, T)

# Genera un mapa de calor para un subconjunto de los datos.
mSet<-PlotSubHeatMap(mSet, "heatmap_2_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, 10.0,0.02,10, 10, "tanova", 25, T, T, T, F, T, T, T,T,"overview")

# Cálculo del fold change, para la identificación de variables con cambios notables entre los grupos.
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)

# Genera un grafico de Fold change.
mSet<-PlotFC(mSet, "fc_1_", "png", 72, width=NA)

# realiza pruebas t para comparar la media de variables entre dos grupos.
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "raw", FALSE)

# Genera un grafico de los resultados de las pruebas t.
mSet<-PlotTT(mSet, "tt_2_", "png", 72, width=NA)


