# SummarizedExperiment.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment")
}

if (!requireNamespace("S4Vectors", quietly = TRUE)) {
  BiocManager::install("S4Vectors")
}
# Se cargan las librerías necesarias

library(ggplot2)
library(S4Vectors)
library(SummarizedExperiment)

### ----------------------------------------------------------
### 1. PREPARACIÓN DE LOS COMPONENTES DEL SummarizedExperiment
### ----------------------------------------------------------

# Se leen los datos del archivo CSV

datos <- read.csv("c:/PEC1/ADO/data/microbiota_SummarizedExperiment.csv")

### a) assays: Los datos principales (matriz de metabolitos)

# Se extrae la matriz de expresión de metabolitos (todas las columnas excepto las dos primeras)

assay_data <- as.matrix(datos[, -c(1, 2)])
rownames(assay_data) <- datos$Muestra  # Las muestras son las filas
colnames(assay_data) <- colnames(datos)[-c(1, 2)]  # Los metabolitos son las columnas menos las 2 primeras
colnames(assay_data)

# En SummarizedExperiment, las muestras deben ser columnas y las características filas
# por lo que se transpone la matriz

assay_data <- t(assay_data)


# Se crea una lista de assays (se pueden tener múltiples matrices de datos)

assays_list <- list(
  raw_counts = assay_data          # Valores originales
)

### b) colData: Metadatos de las muestras (columnas)
# Se crea un DataFrame con información sobre cada muestra


col_data <- DataFrame(
  Muestra = datos$Muestra,
  Grupo = datos$Grupo,
  row.names = datos$Muestra
)


### c) rowData: Metadatos de las características (metabolitos)
# Se crea un DataFrame con información sobre cada metabolito

row_data <- DataFrame(
  Metabolito = colnames(datos)[-c(1, 2)],
  Clase = ifelse(grepl("^[0-9]", colnames(datos)[-c(1, 2)]), 
                 "Químico", "Amino_Acido"),  # Ejemplo de clasificación ficticia
  row.names = colnames(datos)[-c(1, 2)]
)

### d) metadata: Información adicional sobre el experimento
# Se crea una lista con metadatos generales del experimento

experiment_metadata <- list(
  Autor = "Salvador Rizo",
  Fecha = Sys.Date(),
  Description = "Datos de metabolómica"
)

### ----------------------------------------------------------
### 2. CONSTRUCCIÓN DEL OBJETO SummarizedExperiment
### ----------------------------------------------------------

# Se crea el objeto SummarizedExperiment combinando todos los componentes

se <- SummarizedExperiment(
  assays = assays_list,      # Lista de matrices de datos
  colData = col_data,        # Metadatos de las muestras
  rowData = row_data,        # Metadatos de las características
  metadata = experiment_metadata # Metadatos del experimento
)


### ----------------------------------------------------------
### 3. ACCESO A LOS COMPONENTES DEL OBJETO
### ----------------------------------------------------------

### a) Acceso a los assays
# Se obtiene la matriz de datos crudos

raw_data <- assay(se, "raw_counts")


### b) Se accede a los metadatos de las muestras (colData)
muestras_info <- colData(se)

### c) Se accede a los metadatos de las características (rowData)

caracteristicas_info <- rowData(se)

### d) Se accede a los metadatos del experimento (metadata)

experimento_info <- metadata(se)


# Guarda el objeto en un archivo .Rda
save(se, file = "SummarizedExperiment.Rda")
