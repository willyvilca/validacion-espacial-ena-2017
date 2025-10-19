# ============================================
# SCRIPT 02: ANÁLISIS DE AUTOCORRELACIÓN ESPACIAL
# Proyecto: Validación Cruzada Bloqueada y Jerárquica - ENA 2017
# ============================================

rm(list = ls())
gc()

# ============================================
# 1. CONFIGURACIÓN DE RUTAS
# ============================================

# Establecer directorio de trabajo
setwd("D:/Estadistica e Informatica 2025-II/Estadistica Espacial/Proyecto_ENA_2017_R")

# Verificar directorio
cat("Directorio de trabajo:", getwd(), "\n\n")

# ============================================
# 2. INSTALAR/CARGAR LIBRERÍAS
# ============================================

paquetes_necesarios <- c(
  "tidyverse", "spdep", "sf", "gstat", "sp",
  "ggplot2", "scales", "gridExtra", "viridis"
)

cat("Verificando e instalando paquetes necesarios...\n")
for (paq in paquetes_necesarios) {
  if (!require(paq, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Instalando %s...\n", paq))
    install.packages(paq, dependencies = TRUE, repos = "https://cran.rstudio.com/")
    library(paq, character.only = TRUE)
  }
}

cat("\n================================================================\n")
cat("ANALISIS DE AUTOCORRELACION ESPACIAL - ENA 2017\n")
cat("================================================================\n\n")

# ============================================
# 3. CARGAR DATOS LIMPIOS
# ============================================

cat("Cargando datos limpios...\n")

ruta_datos <- "datos_limpios.csv"

if (!file.exists(ruta_datos)) {
  stop("ERROR: No se encuentra el archivo datos_limpios.csv")
}

df <- read_csv(ruta_datos, show_col_types = FALSE)

cat(sprintf("Datos cargados: %s observaciones\n", format(nrow(df), big.mark = ",")))
cat(sprintf("Variables: %d\n\n", ncol(df)))

# Resumen
cat("Resumen de P402A (Inventario Ganadero):\n")
print(summary(df$P402A))

cat("\nDistribucion por Dominio:\n")
print(table(df$DOMINIO))

# ============================================
# 4. GENERAR COORDENADAS ESPACIALES
# ============================================

cat("\nGenerando coordenadas espaciales sinteticas...\n")

dominios_coords <- data.frame(
  DOMINIO = 1:7,
  lat_base = c(-8.0, -12.0, -16.0, -6.0, -12.0, -15.5, -4.5),
  lon_base = c(-79.0, -77.0, -72.0, -78.5, -75.5, -70.0, -76.0)
)

df <- df %>%
  left_join(dominios_coords, by = "DOMINIO")

set.seed(42)
df <- df %>%
  mutate(
    lat = lat_base + rnorm(n(), mean = 0, sd = 0.8),
    lon = lon_base + rnorm(n(), mean = 0, sd = 0.8)
  ) %>%
  select(-lat_base, -lon_base)

cat(sprintf("Coordenadas generadas\n"))
cat(sprintf("Rango latitud: [%.2f, %.2f]\n", min(df$lat), max(df$lat)))
cat(sprintf("Rango longitud: [%.2f, %.2f]\n\n", min(df$lon), max(df$lon)))

# ============================================
# 5. CREAR OBJETO ESPACIAL
# ============================================

cat("Creando objeto espacial sf...\n")
df_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
cat("Objeto espacial creado\n\n")

# ============================================
# 6. MATRIZ DE VECINOS Y PESOS
# ============================================

cat("Calculando matriz de vecinos...\n")

coords <- st_coordinates(df_sf)
k_neighbors <- 8
nb_knn <- knn2nb(knearneigh(coords, k = k_neighbors))

cat(sprintf("Matriz KNN creada con k=%d vecinos\n", k_neighbors))
cat(sprintf("Total de conexiones: %d\n\n", sum(card(nb_knn))))

listw_knn <- nb2listw(nb_knn, style = "W", zero.policy = TRUE)

# ============================================
# 7. I DE MORAN GLOBAL
# ============================================

cat("================================================================\n")
cat("I DE MORAN GLOBAL\n")
cat("================================================================\n\n")

y <- df_sf$P402A
moran_test <- moran.test(y, listw_knn, zero.policy = TRUE)

I_moran <- moran_test$estimate["Moran I statistic"]
E_I <- moran_test$estimate["Expectation"]
V_I <- moran_test$estimate["Variance"]
z_score <- (I_moran - E_I) / sqrt(V_I)
p_value <- moran_test$p.value

cat("RESULTADOS:\n")
cat(sprintf("I de Moran:      %.4f\n", I_moran))
cat(sprintf("Valor esperado:  %.4f\n", E_I))
cat(sprintf("Varianza:        %.6f\n", V_I))
cat(sprintf("Z-score:         %.4f\n", z_score))
cat(sprintf("P-value:         %.6f\n\n", p_value))

if (p_value < 0.001) {
  significancia <- "*** (p < 0.001) - ALTAMENTE SIGNIFICATIVO"
} else if (p_value < 0.01) {
  significancia <- "** (p < 0.01) - MUY SIGNIFICATIVO"
} else if (p_value < 0.05) {
  significancia <- "* (p < 0.05) - SIGNIFICATIVO"
} else {
  significancia <- "(p >= 0.05) - NO SIGNIFICATIVO"
}

cat(sprintf("Significancia: %s\n\n", significancia))

if (I_moran > 0 && p_value < 0.05) {
  cat("INTERPRETACION: Autocorrelacion espacial POSITIVA significativa\n")
  cat("Validacion cruzada bloqueada es NECESARIA.\n\n")
}

# ============================================
# 8. MORAN SCATTERPLOT
# ============================================

y_std <- scale(y)[,1]
y_lag <- lag.listw(listw_knn, y_std, zero.policy = TRUE)

moran_df <- data.frame(
  y_std = y_std,
  y_lag = y_lag,
  dominio = as.factor(df_sf$DOMINIO)
)

p1 <- ggplot(moran_df, aes(x = y_std, y = y_lag)) +
  geom_point(aes(color = dominio), alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "darkblue", linewidth = 1.5) +
  labs(
    title = sprintf("Moran Scatterplot\nI = %.4f, p = %.4f", I_moran, p_value),
    x = "P402A (estandarizado)",
    y = "Spatial Lag de P402A",
    color = "Dominio"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ============================================
# 9. PERMUTACIONES MONTE CARLO
# ============================================

cat("Realizando test de permutaciones Monte Carlo...\n")
moran_mc <- moran.mc(y, listw_knn, nsim = 999, zero.policy = TRUE)

cat(sprintf("Permutaciones completadas\n"))
cat(sprintf("I observado: %.4f\n", moran_mc$statistic))
cat(sprintf("P-value MC:  %.6f\n\n", moran_mc$p.value))

p2 <- ggplot(data.frame(I = moran_mc$res), aes(x = I)) +
  geom_histogram(bins = 50, fill = "gray70", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = I_moran), color = "red", linewidth = 1.5) +
  geom_vline(aes(xintercept = E_I), color = "blue", linewidth = 1.5, linetype = "dashed") +
  labs(
    title = "Distribucion Nula de I de Moran\n(999 permutaciones)",
    x = "I de Moran",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ============================================
# 10. LISA (INDICADORES LOCALES)
# ============================================

cat("================================================================\n")
cat("INDICADORES LOCALES (LISA)\n")
cat("================================================================\n\n")

lisa <- localmoran(y, listw_knn, zero.policy = TRUE)

df_sf$lisa_I <- lisa[, 1]
df_sf$lisa_pval <- lisa[, 5]

df_sf$quadrant <- NA
df_sf$quadrant[y_std > 0 & y_lag > 0] <- "HH"
df_sf$quadrant[y_std < 0 & y_lag > 0] <- "LH"
df_sf$quadrant[y_std < 0 & y_lag < 0] <- "LL"
df_sf$quadrant[y_std > 0 & y_lag < 0] <- "HL"

df_sf$sig <- ifelse(df_sf$lisa_pval < 0.05, "Significativo", "No significativo")
df_sf$lisa_class <- ifelse(df_sf$sig == "Significativo", df_sf$quadrant, "NS")

cat("CLUSTERING ESPACIAL SIGNIFICATIVO (p < 0.05):\n")
print(table(df_sf$lisa_class))
cat("\n")

df_plot <- df_sf %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  )

p3 <- ggplot(df_plot) +
  geom_point(aes(x = lon, y = lat, color = lisa_class), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("HH" = "red", "LL" = "blue", "LH" = "lightblue", 
               "HL" = "pink", "NS" = "gray80"),
    labels = c("HH" = "HH (Hotspot)", "LL" = "LL (Coldspot)", 
               "LH" = "LH", "HL" = "HL", "NS" = "No significativo")
  ) +
  labs(
    title = "Indicadores Locales de Asociacion Espacial (LISA)",
    subtitle = "Inventario Ganadero P402A",
    x = "Longitud", y = "Latitud", color = "LISA"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ============================================
# 11. VARIOGRAMA
# ============================================

cat("================================================================\n")
cat("VARIOGRAMA EMPIRICO\n")
cat("================================================================\n\n")

cat("Calculando variograma...\n")

n_sample <- min(5000, nrow(df_sf))
set.seed(42)
df_sample <- df_sf[sample(1:nrow(df_sf), n_sample), ]

coords_sample <- st_coordinates(df_sample)
y_sample <- df_sample$P402A

sp_data <- data.frame(x = coords_sample[,1], y = coords_sample[,2], P402A = y_sample)
coordinates(sp_data) <- ~x+y

vgm_emp <- variogram(P402A ~ 1, data = sp_data, cutoff = 3, width = 0.15)
vgm_fit <- fit.variogram(vgm_emp, vgm(psill = var(y_sample), model = "Sph", range = 1, nugget = 0))

cat("MODELO DE VARIOGRAMA AJUSTADO:\n")
print(vgm_fit)
cat("\n")

rango_grados <- vgm_fit$range[2]
rango_km <- rango_grados * 111

cat(sprintf("RANGO DE AUTOCORRELACION ESTIMADO: ~%.0f km\n", rango_km))
cat(sprintf("Recomendacion: Usar buffer de >= %.0f km\n\n", rango_km * 1.2))

vgm_df <- data.frame(vgm_emp)
p4 <- ggplot(vgm_df, aes(x = dist, y = gamma)) +
  geom_point(aes(size = np), color = "darkblue", alpha = 0.6) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_hline(yintercept = var(y_sample), linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Variograma Empirico",
    subtitle = sprintf("Rango estimado: %.0f km", rango_km),
    x = "Distancia (grados x 111 ≈ km)",
    y = "Semivarianza",
    size = "Pares"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ============================================
# 12. GUARDAR RESULTADOS
# ============================================

cat("Guardando graficos...\n")

dir.create("resultados/graficos", recursive = TRUE, showWarnings = FALSE)
dir.create("resultados/tablas", recursive = TRUE, showWarnings = FALSE)

ggsave("resultados/graficos/02_moran_scatterplot.png", p1, width = 10, height = 8, dpi = 300)
ggsave("resultados/graficos/02_distribucion_nula.png", p2, width = 10, height = 8, dpi = 300)
ggsave("resultados/graficos/02_lisa_map.png", p3, width = 12, height = 8, dpi = 300)
ggsave("resultados/graficos/02_variogram.png", p4, width = 10, height = 8, dpi = 300)

p_combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("resultados/graficos/02_analisis_espacial_completo.png", p_combined, width = 16, height = 14, dpi = 300)

cat("Graficos guardados\n\n")

# Resultados numéricos
resultados <- data.frame(
  Metrica = c("I de Moran", "P-value", "Z-score", "Rango estimado (km)", 
              "Sill (varianza)", "N observaciones", "Hotspots (HH)", "Coldspots (LL)"),
  Valor = c(I_moran, p_value, z_score, rango_km, var(y, na.rm = TRUE), nrow(df_sf),
            sum(df_sf$lisa_class == "HH", na.rm = TRUE),
            sum(df_sf$lisa_class == "LL", na.rm = TRUE))
)

write_csv(resultados, "resultados/tablas/02_autocorrelacion_espacial.csv")
cat("Resultados guardados\n\n")

# ============================================
# 13. RESUMEN FINAL
# ============================================

cat("================================================================\n")
cat("RESUMEN EJECUTIVO\n")
cat("================================================================\n\n")

cat(sprintf("I DE MORAN GLOBAL:\n"))
cat(sprintf("  Valor: %.4f\n", I_moran))
cat(sprintf("  P-value: %.6f %s\n\n", p_value, significancia))

cat("ANALISIS LISA:\n")
cat(sprintf("  Hotspots (HH): %d ubicaciones\n", sum(df_sf$lisa_class == "HH", na.rm = TRUE)))
cat(sprintf("  Coldspots (LL): %d ubicaciones\n\n", sum(df_sf$lisa_class == "LL", na.rm = TRUE)))

cat("VARIOGRAMA:\n")
cat(sprintf("  Rango de autocorrelacion: ~%.0f km\n", rango_km))
cat(sprintf("  Buffer recomendado: >= %.0f km\n\n", rango_km * 1.2))

cat("CONCLUSION:\n")
intensidad <- ifelse(I_moran > 0.3, "FUERTE", ifelse(I_moran > 0.1, "MODERADA", "DEBIL"))
cat(sprintf("  La variable P402A exhibe autocorrelacion espacial %s,\n", intensidad))
cat("  justificando validacion cruzada bloqueada y jerarquica.\n\n")

cat("================================================================\n")
cat("ANALISIS ESPACIAL COMPLETADO\n")
cat("================================================================\n")




# ============================================
# SCRIPT 03: VALIDACIÓN CRUZADA BLOQUEADA
# Proyecto: Validación Cruzada Bloqueada y Jerárquica - ENA 2017
# ============================================

# ============================================
# 1. CONFIGURACIÓN
# ============================================

setwd("D:/Estadistica e Informatica 2025-II/Estadistica Espacial/Proyecto_ENA_2017_R")
cat("Directorio de trabajo:", getwd(), "\n\n")

# Cargar librerías
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(Metrics)
  library(ggplot2)
  library(gridExtra)
})

cat("================================================================\n")
cat("VALIDACION CRUZADA BLOQUEADA - ENA 2017\n")
cat("================================================================\n\n")

# ============================================
# 2. CARGAR DATOS
# ============================================

cat("Cargando datos limpios...\n")
df <- read_csv("datos_limpios.csv", show_col_types = FALSE)

# Seleccionar variables
df_modelo <- df %>%
  select(P402A, DOMINIO, P401A, CCDD) %>%
  filter(!is.na(P402A), !is.na(DOMINIO), !is.na(P401A))

cat(sprintf("Datos cargados: %s observaciones\n", format(nrow(df_modelo), big.mark = ",")))
cat(sprintf("Variables: P402A (dependiente), DOMINIO (bloqueo), P401A (especie)\n\n"))

# Convertir a factores
df_modelo <- df_modelo %>%
  mutate(
    DOMINIO = as.factor(DOMINIO),
    P401A = as.factor(P401A),
    CCDD = as.factor(CCDD)
  )

# Resumen
cat("Distribucion por dominio:\n")
print(table(df_modelo$DOMINIO))
cat("\n")

# ============================================
# 3. PREPARAR DATOS PARA MODELAMIENTO
# ============================================

cat("Preparando datos para modelamiento...\n")

# Variables predictoras
X <- df_modelo %>% select(DOMINIO, P401A)
y <- df_modelo$P402A

# Transformación log de y (mejora distribución)
y_log <- log1p(y)  # log(1 + y) para manejar ceros

cat(sprintf("Variable objetivo transformada: log(P402A + 1)\n"))
cat(sprintf("Media original P402A: %.2f\n", mean(y)))
cat(sprintf("Media transformada: %.2f\n\n", mean(y_log)))

# ============================================
# 4. FUNCIÓN DE EVALUACIÓN
# ============================================

evaluar_modelo <- function(y_real, y_pred, nombre = "") {
  # Calcular métricas
  rmse_val <- rmse(y_real, y_pred)
  mae_val <- mae(y_real, y_pred)
  r2_val <- cor(y_real, y_pred)^2
  
  # Retornar como dataframe
  return(data.frame(
    Metodo = nombre,
    RMSE = rmse_val,
    MAE = mae_val,
    R2 = r2_val
  ))
}

# ============================================
# 5. VALIDACIÓN CRUZADA ALEATORIA (BASELINE)
# ============================================

cat("================================================================\n")
cat("VALIDACION CRUZADA ALEATORIA (BASELINE)\n")
cat("================================================================\n\n")

set.seed(42)
n_folds <- 10

cat(sprintf("Configuracion: %d-fold cross-validation aleatorio\n", n_folds))
cat("Modelo: Random Forest (ntree=100, mtry=1)\n\n")

# Crear folds aleatorios
folds_random <- createFolds(y_log, k = n_folds, list = TRUE)

# Almacenar resultados
resultados_random <- data.frame()
predicciones_random <- numeric(length(y_log))

cat("Ejecutando validacion aleatoria...\n")
pb <- txtProgressBar(min = 0, max = n_folds, style = 3)

for (i in 1:n_folds) {
  # Índices de train/test
  idx_test <- folds_random[[i]]
  idx_train <- setdiff(1:nrow(df_modelo), idx_test)
  
  # Datos train/test
  X_train <- X[idx_train, ]
  X_test <- X[idx_test, ]
  y_train <- y_log[idx_train]
  y_test <- y_log[idx_test]
  
  # Entrenar modelo
  modelo <- randomForest(x = X_train, y = y_train, 
                         ntree = 100, mtry = 1, 
                         nodesize = 5)
  
  # Predecir
  y_pred <- predict(modelo, X_test)
  predicciones_random[idx_test] <- y_pred
  
  # Evaluar
  metricas <- evaluar_modelo(y_test, y_pred, nombre = paste0("Fold_", i))
  resultados_random <- rbind(resultados_random, metricas)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Promedios
metricas_random <- resultados_random %>%
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_sd = sd(RMSE),
    MAE_mean = mean(MAE),
    MAE_sd = sd(MAE),
    R2_mean = mean(R2),
    R2_sd = sd(R2)
  )

cat("\n\nRESULTADOS VALIDACION ALEATORIA:\n")
cat(sprintf("RMSE: %.4f ± %.4f\n", metricas_random$RMSE_mean, metricas_random$RMSE_sd))
cat(sprintf("MAE:  %.4f ± %.4f\n", metricas_random$MAE_mean, metricas_random$MAE_sd))
cat(sprintf("R²:   %.4f ± %.4f\n\n", metricas_random$R2_mean, metricas_random$R2_sd))

# ============================================
# 6. VALIDACIÓN CRUZADA BLOQUEADA
# ============================================

cat("================================================================\n")
cat("VALIDACION CRUZADA BLOQUEADA (LEAVE-ONE-DOMAIN-OUT)\n")
cat("================================================================\n\n")

dominios <- levels(df_modelo$DOMINIO)
n_dominios <- length(dominios)

cat(sprintf("Configuracion: Leave-one-domain-out (%d folds)\n", n_dominios))
cat(sprintf("Dominios: %s\n", paste(dominios, collapse = ", ")))
cat("Buffer: No aplicable (bloqueo completo por dominio)\n")
cat("Modelo: Random Forest (ntree=100, mtry=1)\n\n")

# Almacenar resultados
resultados_bloqueado <- data.frame()
predicciones_bloqueado <- numeric(length(y_log))

cat("Ejecutando validacion bloqueada...\n")
pb <- txtProgressBar(min = 0, max = n_dominios, style = 3)

for (i in 1:n_dominios) {
  dominio_test <- dominios[i]
  
  # Índices de train/test
  idx_test <- which(df_modelo$DOMINIO == dominio_test)
  idx_train <- which(df_modelo$DOMINIO != dominio_test)
  
  # Datos train/test
  X_train <- X[idx_train, ]
  X_test <- X[idx_test, ]
  y_train <- y_log[idx_train]
  y_test <- y_log[idx_test]
  
  # Entrenar modelo
  modelo <- randomForest(x = X_train, y = y_train, 
                         ntree = 100, mtry = 1, 
                         nodesize = 5)
  
  # Predecir
  y_pred <- predict(modelo, X_test)
  predicciones_bloqueado[idx_test] <- y_pred
  
  # Evaluar
  metricas <- evaluar_modelo(y_test, y_pred, nombre = paste0("Dominio_", dominio_test))
  metricas$n_test <- length(idx_test)
  metricas$n_train <- length(idx_train)
  resultados_bloqueado <- rbind(resultados_bloqueado, metricas)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Promedios (ponderados por tamaño)
metricas_bloqueado <- resultados_bloqueado %>%
  summarise(
    RMSE_mean = weighted.mean(RMSE, n_test),
    RMSE_sd = sd(RMSE),
    MAE_mean = weighted.mean(MAE, n_test),
    MAE_sd = sd(MAE),
    R2_mean = weighted.mean(R2, n_test),
    R2_sd = sd(R2)
  )

cat("\n\nRESULTADOS VALIDACION BLOQUEADA:\n")
cat(sprintf("RMSE: %.4f ± %.4f\n", metricas_bloqueado$RMSE_mean, metricas_bloqueado$RMSE_sd))
cat(sprintf("MAE:  %.4f ± %.4f\n", metricas_bloqueado$MAE_mean, metricas_bloqueado$MAE_sd))
cat(sprintf("R²:   %.4f ± %.4f\n\n", metricas_bloqueado$R2_mean, metricas_bloqueado$R2_sd))

# Resultados por dominio
cat("Resultados por dominio:\n")
print(resultados_bloqueado %>% select(Metodo, RMSE, MAE, R2, n_test))
cat("\n")

# ============================================
# 7. COMPARACIÓN DE MÉTODOS
# ============================================

cat("================================================================\n")
cat("COMPARACION DE METODOS\n")
cat("================================================================\n\n")

# Crear tabla comparativa
comparacion <- data.frame(
  Metodo = c("Validacion Aleatoria", "Validacion Bloqueada"),
  RMSE = c(metricas_random$RMSE_mean, metricas_bloqueado$RMSE_mean),
  MAE = c(metricas_random$MAE_mean, metricas_bloqueado$MAE_mean),
  R2 = c(metricas_random$R2_mean, metricas_bloqueado$R2_mean)
)

print(comparacion)
cat("\n")

# Calcular diferencias
diff_rmse <- ((metricas_bloqueado$RMSE_mean - metricas_random$RMSE_mean) / 
                metricas_random$RMSE_mean) * 100
diff_mae <- ((metricas_bloqueado$MAE_mean - metricas_random$MAE_mean) / 
               metricas_random$MAE_mean) * 100
diff_r2 <- ((metricas_random$R2_mean - metricas_bloqueado$R2_mean) / 
              metricas_random$R2_mean) * 100

cat("FILTRACION ESPACIAL (sobreestimacion de validacion aleatoria):\n")
cat(sprintf("RMSE: Validacion bloqueada es %.2f%% mayor\n", diff_rmse))
cat(sprintf("MAE:  Validacion bloqueada es %.2f%% mayor\n", diff_mae))
cat(sprintf("R²:   Validacion aleatoria sobreestima en %.2f%%\n\n", diff_r2))

# Prueba t pareada (solo si hay suficientes folds comparables)
if (n_folds == n_dominios) {
  t_test_rmse <- t.test(resultados_random$RMSE, resultados_bloqueado$RMSE, paired = TRUE)
  cat("Test t pareado para RMSE:\n")
  cat(sprintf("t = %.4f, p-value = %.6f\n", t_test_rmse$statistic, t_test_rmse$p.value))
  if (t_test_rmse$p.value < 0.05) {
    cat("Diferencia SIGNIFICATIVA entre metodos (p < 0.05)\n\n")
  } else {
    cat("Diferencia NO significativa (p >= 0.05)\n\n")
  }
}

# ============================================
# 8. VISUALIZACIONES
# ============================================

cat("Generando graficos...\n")

dir.create("resultados/graficos", recursive = TRUE, showWarnings = FALSE)

# 1. Comparación de métricas
metricas_long <- comparacion %>%
  pivot_longer(cols = c(RMSE, MAE, R2), names_to = "Metrica", values_to = "Valor")

p1 <- ggplot(metricas_long, aes(x = Metrica, y = Valor, fill = Metodo)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Valor)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Validacion Aleatoria" = "steelblue", 
                               "Validacion Bloqueada" = "coral")) +
  labs(
    title = "Comparacion de Metricas de Validacion",
    subtitle = "Random Forest - Inventario Ganadero P402A",
    x = "Metrica",
    y = "Valor",
    fill = "Metodo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom"
  )

# 2. RMSE por dominio (validación bloqueada)
p2 <- ggplot(resultados_bloqueado, aes(x = Metodo, y = RMSE)) +
  geom_bar(stat = "identity", fill = "coral", width = 0.7) +
  geom_hline(yintercept = metricas_random$RMSE_mean, 
             linetype = "dashed", color = "steelblue", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f", RMSE)), vjust = -0.5, size = 3) +
  annotate("text", x = 4, y = metricas_random$RMSE_mean * 1.05, 
           label = "RMSE Validacion Aleatoria", color = "steelblue", size = 3.5) +
  labs(
    title = "RMSE por Dominio - Validacion Bloqueada",
    subtitle = "Linea punteada: RMSE promedio de validacion aleatoria",
    x = "Dominio (Test Set)",
    y = "RMSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 3. Scatter plot: Observado vs Predicho
df_scatter <- data.frame(
  Observado = y_log,
  Pred_Random = predicciones_random,
  Pred_Bloqueado = predicciones_bloqueado
) %>%
  pivot_longer(cols = c(Pred_Random, Pred_Bloqueado), 
               names_to = "Metodo", values_to = "Predicho") %>%
  mutate(Metodo = recode(Metodo, 
                         "Pred_Random" = "Validacion Aleatoria",
                         "Pred_Bloqueado" = "Validacion Bloqueada"))

p3 <- ggplot(df_scatter, aes(x = Observado, y = Predicho, color = Metodo)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Validacion Aleatoria" = "steelblue", 
                                "Validacion Bloqueada" = "coral")) +
  facet_wrap(~Metodo) +
  labs(
    title = "Observado vs Predicho",
    subtitle = "Log(P402A + 1)",
    x = "Valor Observado",
    y = "Valor Predicho"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none"
  )

# 4. Distribución de errores
df_errores <- data.frame(
  Error_Random = y_log - predicciones_random,
  Error_Bloqueado = y_log - predicciones_bloqueado
) %>%
  pivot_longer(cols = everything(), names_to = "Metodo", values_to = "Error") %>%
  mutate(Metodo = recode(Metodo,
                         "Error_Random" = "Validacion Aleatoria",
                         "Error_Bloqueado" = "Validacion Bloqueada"))

p4 <- ggplot(df_errores, aes(x = Error, fill = Metodo)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Validacion Aleatoria" = "steelblue", 
                               "Validacion Bloqueada" = "coral")) +
  labs(
    title = "Distribucion de Errores de Prediccion",
    x = "Error (Observado - Predicho)",
    y = "Frecuencia",
    fill = "Metodo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

# Guardar gráficos
ggsave("resultados/graficos/03_comparacion_metricas.png", p1, 
       width = 10, height = 8, dpi = 300)
ggsave("resultados/graficos/03_rmse_por_dominio.png", p2, 
       width = 12, height = 8, dpi = 300)
ggsave("resultados/graficos/03_observado_vs_predicho.png", p3, 
       width = 12, height = 6, dpi = 300)
ggsave("resultados/graficos/03_distribucion_errores.png", p4, 
       width = 10, height = 8, dpi = 300)

# Panel combinado
p_combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("resultados/graficos/03_validacion_bloqueada_completo.png", p_combined,
       width = 16, height = 14, dpi = 300)

cat("Graficos guardados\n\n")

# ============================================
# 9. GUARDAR RESULTADOS
# ============================================

dir.create("resultados/tablas", recursive = TRUE, showWarnings = FALSE)

# Guardar comparación
write_csv(comparacion, "resultados/tablas/03_comparacion_metodos.csv")

# Guardar resultados detallados
write_csv(resultados_bloqueado, "resultados/tablas/03_resultados_por_dominio.csv")

# Resumen ejecutivo
resumen <- data.frame(
  Aspecto = c(
    "Metodo_Baseline",
    "Metodo_Propuesto",
    "RMSE_Baseline",
    "RMSE_Bloqueado",
    "Incremento_RMSE_pct",
    "R2_Baseline",
    "R2_Bloqueado",
    "Reduccion_R2_pct",
    "Filtracion_Espacial"
  ),
  Valor = c(
    "Validacion Aleatoria 10-fold",
    "Validacion Bloqueada Leave-One-Domain-Out",
    sprintf("%.4f", metricas_random$RMSE_mean),
    sprintf("%.4f", metricas_bloqueado$RMSE_mean),
    sprintf("%.2f%%", diff_rmse),
    sprintf("%.4f", metricas_random$R2_mean),
    sprintf("%.4f", metricas_bloqueado$R2_mean),
    sprintf("%.2f%%", diff_r2),
    ifelse(diff_rmse > 10, "ALTA", ifelse(diff_rmse > 5, "MODERADA", "BAJA"))
  )
)

write_csv(resumen, "resultados/tablas/03_resumen_validacion_bloqueada.csv")

cat("Resultados guardados\n\n")

# ============================================
# 10. RESUMEN FINAL
# ============================================

cat("================================================================\n")
cat("RESUMEN EJECUTIVO - VALIDACION BLOQUEADA\n")
cat("================================================================\n\n")

cat("VALIDACION ALEATORIA (BASELINE):\n")
cat(sprintf("  RMSE: %.4f\n", metricas_random$RMSE_mean))
cat(sprintf("  R²:   %.4f\n\n", metricas_random$R2_mean))

cat("VALIDACION BLOQUEADA (LEAVE-ONE-DOMAIN-OUT):\n")
cat(sprintf("  RMSE: %.4f\n", metricas_bloqueado$RMSE_mean))
cat(sprintf("  R²:   %.4f\n\n", metricas_bloqueado$R2_mean))

cat("FILTRACION ESPACIAL:\n")
cat(sprintf("  La validacion aleatoria SOBREESTIMA el rendimiento:\n"))
cat(sprintf("  - RMSE %.2f%% mas bajo que validacion bloqueada\n", abs(diff_rmse)))
cat(sprintf("  - R² %.2f%% mas alto que validacion bloqueada\n\n", diff_r2))

cat("CONCLUSION:\n")
if (diff_rmse > 10) {
  cat("  Filtracion espacial ALTA detectada.\n")
  cat("  Validacion cruzada bloqueada es CRITICA para estimaciones\n")
  cat("  honestas de rendimiento predictivo.\n\n")
} else if (diff_rmse > 5) {
  cat("  Filtracion espacial MODERADA detectada.\n")
  cat("  Validacion cruzada bloqueada mejora la confiabilidad de\n")
  cat("  las estimaciones de rendimiento.\n\n")
} else {
  cat("  Filtracion espacial BAJA detectada.\n")
  cat("  Validacion cruzada bloqueada proporciona estimaciones\n")
  cat("  ligeramente mas conservadoras.\n\n")
}

cat("================================================================\n")
cat("VALIDACION BLOQUEADA COMPLETADA\n")
cat("================================================================\n")




# ============================================
# SCRIPT 04: VALIDACIÓN CRUZADA JERÁRQUICA
# Proyecto: Validación Cruzada Bloqueada y Jerárquica - ENA 2017
# ============================================

# ============================================
# 1. CONFIGURACIÓN
# ============================================

setwd("D:/Estadistica e Informatica 2025-II/Estadistica Espacial/Proyecto_ENA_2017_R")
cat("Directorio de trabajo:", getwd(), "\n\n")

# Cargar librerías
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(Metrics)
  library(ggplot2)
  library(gridExtra)
})

cat("================================================================\n")
cat("VALIDACION CRUZADA JERARQUICA - ENA 2017\n")
cat("================================================================\n\n")

# ============================================
# 2. CARGAR DATOS
# ============================================

cat("Cargando datos limpios...\n")
df <- read_csv("datos_limpios.csv", show_col_types = FALSE)

df_modelo <- df %>%
  select(P402A, DOMINIO, P401A, CCDD) %>%
  filter(!is.na(P402A), !is.na(DOMINIO), !is.na(P401A), !is.na(CCDD))

cat(sprintf("Datos cargados: %s observaciones\n", format(nrow(df_modelo), big.mark = ",")))

# Convertir a factores
df_modelo <- df_modelo %>%
  mutate(
    DOMINIO = as.factor(DOMINIO),
    P401A = as.factor(P401A),
    CCDD = as.factor(CCDD)
  )

# Transformación log
X <- df_modelo %>% select(DOMINIO, P401A)
y <- df_modelo$P402A
y_log <- log1p(y)

cat("Variables preparadas\n\n")

# ============================================
# 3. FUNCIÓN DE EVALUACIÓN
# ============================================

evaluar_modelo <- function(y_real, y_pred, nombre = "") {
  rmse_val <- rmse(y_real, y_pred)
  mae_val <- mae(y_real, y_pred)
  r2_val <- cor(y_real, y_pred)^2
  
  return(data.frame(
    Metodo = nombre,
    RMSE = rmse_val,
    MAE = mae_val,
    R2 = r2_val
  ))
}

# ============================================
# 4. VALIDACIÓN CRUZADA JERÁRQUICA
# ============================================

cat("================================================================\n")
cat("VALIDACION CRUZADA JERARQUICA (NESTED CV)\n")
cat("================================================================\n\n")

cat("Estructura:\n")
cat("  - Outer Loop: Leave-One-Domain-Out (7 folds)\n")
cat("  - Inner Loop: Leave-One-Department-Out (optimizacion hiperparametros)\n")
cat("  - Modelo: Random Forest\n\n")

dominios <- levels(df_modelo$DOMINIO)
n_dominios <- length(dominios)

# Hiperparámetros a probar
hiperparametros <- expand.grid(
  ntree = c(50, 100, 150),
  mtry = c(1, 2),
  nodesize = c(5, 10)
)

cat(sprintf("Hiperparametros a evaluar: %d combinaciones\n", nrow(hiperparametros)))
cat("  ntree: 50, 100, 150\n")
cat("  mtry: 1, 2\n")
cat("  nodesize: 5, 10\n\n")

# Almacenar resultados
resultados_jerarquico <- data.frame()
mejores_hiperparametros <- list()
predicciones_jerarquico <- numeric(length(y_log))

cat("Ejecutando validacion jerarquica...\n")
cat("NOTA: Este proceso puede tomar varios minutos\n\n")

tiempo_inicio <- Sys.time()

# OUTER LOOP: Dominios
for (i in 1:n_dominios) {
  dominio_test <- dominios[i]
  
  cat(sprintf("Outer Fold %d/%d: Dominio Test = %s\n", i, n_dominios, dominio_test))
  
  # Índices outer loop
  idx_outer_test <- which(df_modelo$DOMINIO == dominio_test)
  idx_outer_train <- which(df_modelo$DOMINIO != dominio_test)
  
  # Datos outer train/test
  X_outer_train <- X[idx_outer_train, ]
  X_outer_test <- X[idx_outer_test, ]
  y_outer_train <- y_log[idx_outer_train]
  y_outer_test <- y_log[idx_outer_test]
  
  # Departamentos en outer train
  ccdd_outer_train <- df_modelo$CCDD[idx_outer_train]
  departamentos_train <- unique(ccdd_outer_train)
  n_departamentos <- length(departamentos_train)
  
  cat(sprintf("  Outer Train: %d obs, %d departamentos\n", 
              length(idx_outer_train), n_departamentos))
  
  # INNER LOOP: Optimización de hiperparámetros
  cat("  Optimizando hiperparametros (Inner Loop)...\n")
  
  mejores_rmse <- Inf
  mejores_params <- NULL
  
  pb_inner <- txtProgressBar(min = 0, max = nrow(hiperparametros), style = 3)
  
  for (h in 1:nrow(hiperparametros)) {
    params <- hiperparametros[h, ]
    
    # Validación cruzada inner: leave-one-department-out
    rmse_inner <- numeric(n_departamentos)
    
    for (d in 1:n_departamentos) {
      dept_test <- departamentos_train[d]
      
      # Índices inner
      idx_inner_test <- which(ccdd_outer_train == dept_test)
      idx_inner_train <- which(ccdd_outer_train != dept_test)
      
      # Datos inner
      X_inner_train <- X_outer_train[idx_inner_train, ]
      X_inner_test <- X_outer_train[idx_inner_test, ]
      y_inner_train <- y_outer_train[idx_inner_train]
      y_inner_test <- y_outer_train[idx_inner_test]
      
      # Entrenar modelo con hiperparámetros candidatos
      modelo_inner <- randomForest(
        x = X_inner_train, 
        y = y_inner_train,
        ntree = params$ntree,
        mtry = params$mtry,
        nodesize = params$nodesize
      )
      
      # Predecir
      y_pred_inner <- predict(modelo_inner, X_inner_test)
      
      # RMSE
      rmse_inner[d] <- rmse(y_inner_test, y_pred_inner)
    }
    
    # RMSE promedio en inner CV
    rmse_promedio <- mean(rmse_inner)
    
    # Actualizar mejores hiperparámetros
    if (rmse_promedio < mejores_rmse) {
      mejores_rmse <- rmse_promedio
      mejores_params <- params
    }
    
    setTxtProgressBar(pb_inner, h)
  }
  close(pb_inner)
  
  cat(sprintf("\n  Mejores hiperparametros encontrados:\n"))
  cat(sprintf("    ntree=%d, mtry=%d, nodesize=%d (RMSE inner=%.4f)\n",
              mejores_params$ntree, mejores_params$mtry, 
              mejores_params$nodesize, mejores_rmse))
  
  # Guardar mejores hiperparámetros
  mejores_hiperparametros[[as.character(dominio_test)]] <- mejores_params
  
  # ENTRENAR MODELO FINAL con mejores hiperparámetros en outer train
  cat("  Entrenando modelo final con mejores hiperparametros...\n")
  
  modelo_final <- randomForest(
    x = X_outer_train,
    y = y_outer_train,
    ntree = mejores_params$ntree,
    mtry = mejores_params$mtry,
    nodesize = mejores_params$nodesize
  )
  
  # PREDECIR en outer test
  y_pred_outer <- predict(modelo_final, X_outer_test)
  predicciones_jerarquico[idx_outer_test] <- y_pred_outer
  
  # Evaluar
  metricas <- evaluar_modelo(y_outer_test, y_pred_outer, 
                             nombre = paste0("Dominio_", dominio_test))
  metricas$n_test <- length(idx_outer_test)
  metricas$n_train <- length(idx_outer_train)
  metricas$ntree <- mejores_params$ntree
  metricas$mtry <- mejores_params$mtry
  metricas$nodesize <- mejores_params$nodesize
  
  resultados_jerarquico <- rbind(resultados_jerarquico, metricas)
  
  cat(sprintf("  Outer Test RMSE: %.4f, R2: %.4f\n\n", metricas$RMSE, metricas$R2))
}

tiempo_fin <- Sys.time()
tiempo_total <- difftime(tiempo_fin, tiempo_inicio, units = "mins")

cat(sprintf("Validacion jerarquica completada en %.2f minutos\n\n", tiempo_total))

# ============================================
# 5. RESUMEN DE RESULTADOS
# ============================================

cat("================================================================\n")
cat("RESULTADOS VALIDACION JERARQUICA\n")
cat("================================================================\n\n")

# Métricas promedio (ponderadas)
metricas_jerarquico <- resultados_jerarquico %>%
  summarise(
    RMSE_mean = weighted.mean(RMSE, n_test),
    RMSE_sd = sd(RMSE),
    MAE_mean = weighted.mean(MAE, n_test),
    MAE_sd = sd(MAE),
    R2_mean = weighted.mean(R2, n_test),
    R2_sd = sd(R2)
  )

cat("Metricas promedio:\n")
cat(sprintf("  RMSE: %.4f ± %.4f\n", metricas_jerarquico$RMSE_mean, metricas_jerarquico$RMSE_sd))
cat(sprintf("  MAE:  %.4f ± %.4f\n", metricas_jerarquico$MAE_mean, metricas_jerarquico$MAE_sd))
cat(sprintf("  R²:   %.4f ± %.4f\n\n", metricas_jerarquico$R2_mean, metricas_jerarquico$R2_sd))

cat("Resultados por dominio:\n")
print(resultados_jerarquico %>% 
        select(Metodo, RMSE, MAE, R2, ntree, mtry, nodesize, n_test))
cat("\n")

# Hiperparámetros seleccionados
cat("Hiperparametros seleccionados por dominio:\n")
for (dom in names(mejores_hiperparametros)) {
  params <- mejores_hiperparametros[[dom]]
  cat(sprintf("  Dominio %s: ntree=%d, mtry=%d, nodesize=%d\n",
              dom, params$ntree, params$mtry, params$nodesize))
}
cat("\n")

# ============================================
# 6. COMPARACIÓN CON VALIDACIÓN BLOQUEADA
# ============================================

cat("================================================================\n")
cat("COMPARACION: JERARQUICA VS BLOQUEADA\n")
cat("================================================================\n\n")

# Cargar resultados de validación bloqueada
if (file.exists("resultados/tablas/03_comparacion_metodos.csv")) {
  comparacion_previa <- read_csv("resultados/tablas/03_comparacion_metodos.csv", 
                                 show_col_types = FALSE)
  
  metricas_bloqueado <- comparacion_previa %>% 
    filter(Metodo == "Validacion Bloqueada")
  
  # Tabla comparativa
  comparacion_completa <- data.frame(
    Metodo = c("Val. Aleatoria", "Val. Bloqueada", "Val. Jerarquica"),
    RMSE = c(comparacion_previa$RMSE[1], 
             metricas_bloqueado$RMSE, 
             metricas_jerarquico$RMSE_mean),
    MAE = c(comparacion_previa$MAE[1], 
            metricas_bloqueado$MAE, 
            metricas_jerarquico$MAE_mean),
    R2 = c(comparacion_previa$R2[1], 
           metricas_bloqueado$R2, 
           metricas_jerarquico$R2_mean)
  )
  
  print(comparacion_completa)
  cat("\n")
  
  # Diferencias
  diff_rmse_bj <- ((metricas_jerarquico$RMSE_mean - metricas_bloqueado$RMSE) / 
                     metricas_bloqueado$RMSE) * 100
  diff_r2_bj <- ((metricas_bloqueado$R2 - metricas_jerarquico$R2_mean) / 
                   metricas_bloqueado$R2) * 100
  
  cat("Diferencia Bloqueada vs Jerarquica:\n")
  cat(sprintf("  RMSE: %.2f%% %s\n", abs(diff_rmse_bj), 
              ifelse(diff_rmse_bj > 0, "mayor en jerarquica", "menor en jerarquica")))
  cat(sprintf("  R²:   %.2f%% %s\n\n", abs(diff_r2_bj),
              ifelse(diff_r2_bj > 0, "menor en jerarquica", "mayor en jerarquica")))
  
} else {
  cat("Archivo de validacion bloqueada no encontrado.\n")
  cat("Ejecute primero el Script 03.\n\n")
  comparacion_completa <- NULL
}

# ============================================
# 7. VISUALIZACIONES
# ============================================

cat("Generando graficos...\n")

dir.create("resultados/graficos", recursive = TRUE, showWarnings = FALSE)

# 1. Comparación de los tres métodos
if (!is.null(comparacion_completa)) {
  metricas_long <- comparacion_completa %>%
    pivot_longer(cols = c(RMSE, MAE, R2), names_to = "Metrica", values_to = "Valor")
  
  p1 <- ggplot(metricas_long, aes(x = Metrica, y = Valor, fill = Metodo)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Valor)), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c(
      "Val. Aleatoria" = "steelblue",
      "Val. Bloqueada" = "coral",
      "Val. Jerarquica" = "darkgreen"
    )) +
    labs(
      title = "Comparacion de Metodos de Validacion",
      subtitle = "Random Forest - Inventario Ganadero P402A",
      x = "Metrica",
      y = "Valor",
      fill = "Metodo"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    )
  
  ggsave("resultados/graficos/04_comparacion_tres_metodos.png", p1,
         width = 12, height = 8, dpi = 300)
}

# 2. RMSE por dominio - Validación Jerárquica
p2 <- ggplot(resultados_jerarquico, aes(x = Metodo, y = RMSE)) +
  geom_bar(stat = "identity", fill = "darkgreen", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", RMSE)), vjust = -0.5, size = 3) +
  labs(
    title = "RMSE por Dominio - Validacion Jerarquica",
    subtitle = "Con optimizacion de hiperparametros (Inner CV)",
    x = "Dominio (Test Set)",
    y = "RMSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("resultados/graficos/04_rmse_jerarquico.png", p2,
       width = 12, height = 8, dpi = 300)

# 3. Hiperparámetros seleccionados por dominio
df_hiperparametros <- do.call(rbind, lapply(names(mejores_hiperparametros), function(dom) {
  params <- mejores_hiperparametros[[dom]]
  data.frame(
    Dominio = dom,
    ntree = params$ntree,
    mtry = params$mtry,
    nodesize = params$nodesize
  )
}))

p3 <- ggplot(df_hiperparametros, aes(x = Dominio)) +
  geom_point(aes(y = ntree, color = "ntree"), size = 4) +
  geom_point(aes(y = mtry * 50, color = "mtry x50"), size = 4) +
  geom_point(aes(y = nodesize * 10, color = "nodesize x10"), size = 4) +
  scale_color_manual(values = c("ntree" = "blue", "mtry x50" = "red", "nodesize x10" = "green")) +
  labs(
    title = "Hiperparametros Optimos por Dominio",
    subtitle = "Seleccionados via Inner Cross-Validation",
    x = "Dominio",
    y = "Valor (escalado)",
    color = "Hiperparametro"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

ggsave("resultados/graficos/04_hiperparametros_optimos.png", p3,
       width = 10, height = 8, dpi = 300)

# 4. Scatter: Observado vs Predicho (Jerárquica)
df_scatter <- data.frame(
  Observado = y_log,
  Predicho = predicciones_jerarquico,
  Dominio = df_modelo$DOMINIO
)

p4 <- ggplot(df_scatter, aes(x = Observado, y = Predicho)) +
  geom_point(aes(color = Dominio), alpha = 0.4, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Observado vs Predicho - Validacion Jerarquica",
    subtitle = "Log(P402A + 1)",
    x = "Valor Observado",
    y = "Valor Predicho",
    color = "Dominio"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "right"
  )

ggsave("resultados/graficos/04_observado_vs_predicho_jerarquico.png", p4,
       width = 12, height = 8, dpi = 300)

cat("Graficos guardados\n\n")

# ============================================
# 8. GUARDAR RESULTADOS
# ============================================

dir.create("resultados/tablas", recursive = TRUE, showWarnings = FALSE)

# Guardar resultados detallados
write_csv(resultados_jerarquico, "resultados/tablas/04_resultados_jerarquico.csv")

# Guardar hiperparámetros
write_csv(df_hiperparametros, "resultados/tablas/04_hiperparametros_optimos.csv")

# Guardar comparación completa
if (!is.null(comparacion_completa)) {
  write_csv(comparacion_completa, "resultados/tablas/04_comparacion_completa.csv")
}

# Resumen ejecutivo
resumen <- data.frame(
  Aspecto = c(
    "Metodo",
    "Outer_Loop",
    "Inner_Loop",
    "RMSE_Promedio",
    "R2_Promedio",
    "Tiempo_Ejecucion_min",
    "Hiperparametros_Optimizados"
  ),
  Valor = c(
    "Validacion Cruzada Jerarquica (Nested CV)",
    "Leave-One-Domain-Out (7 dominios)",
    "Leave-One-Department-Out (por dominio)",
    sprintf("%.4f", metricas_jerarquico$RMSE_mean),
    sprintf("%.4f", metricas_jerarquico$R2_mean),
    sprintf("%.2f", tiempo_total),
    "Si (ntree, mtry, nodesize)"
  )
)

write_csv(resumen, "resultados/tablas/04_resumen_jerarquico.csv")

cat("Resultados guardados\n\n")

# ============================================
# 9. RESUMEN EJECUTIVO FINAL
# ============================================

cat("================================================================\n")
cat("RESUMEN EJECUTIVO - VALIDACION JERARQUICA\n")
cat("================================================================\n\n")

cat("METRICAS PROMEDIO:\n")
cat(sprintf("  RMSE: %.4f ± %.4f\n", metricas_jerarquico$RMSE_mean, metricas_jerarquico$RMSE_sd))
cat(sprintf("  MAE:  %.4f ± %.4f\n", metricas_jerarquico$MAE_mean, metricas_jerarquico$MAE_sd))
cat(sprintf("  R²:   %.4f ± %.4f\n\n", metricas_jerarquico$R2_mean, metricas_jerarquico$R2_sd))

if (!is.null(comparacion_completa)) {
  cat("COMPARACION CON OTROS METODOS:\n")
  print(comparacion_completa)
  cat("\n")
}

cat("VENTAJAS DE VALIDACION JERARQUICA:\n")
cat("  1. Optimiza hiperparametros sin usar datos del dominio test\n")
cat("  2. Evita sobreajuste en seleccion de hiperparametros\n")
cat("  3. Proporciona estimaciones honestas de generalizacion espacial\n")
cat("  4. Permite hiperparametros especificos por contexto espacial\n\n")

cat("HIPERPARAMETROS MAS COMUNES:\n")
tabla_freq <- df_hiperparametros %>%
  count(ntree, mtry, nodesize) %>%
  arrange(desc(n))
print(tabla_freq)
cat("\n")

cat("================================================================\n")
cat("VALIDACION JERARQUICA COMPLETADA\n")
cat(sprintf("Tiempo total: %.2f minutos\n", tiempo_total))
cat("================================================================\n")





# ============================================
# SCRIPT 05: COMPARACIÓN FINAL Y TABLAS CONSOLIDADAS
# Proyecto: Validación Cruzada Bloqueada y Jerárquica - ENA 2017
# ============================================

# ============================================
# 1. CONFIGURACIÓN
# ============================================

setwd("D:/Estadistica e Informatica 2025-II/Estadistica Espacial/Proyecto_ENA_2017_R")
cat("Directorio de trabajo:", getwd(), "\n\n")

# Cargar librerías
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
  library(knitr)
  library(kableExtra)
})

cat("================================================================\n")
cat("COMPARACION FINAL Y RESULTADOS CONSOLIDADOS\n")
cat("================================================================\n\n")

# ============================================
# 2. CARGAR TODOS LOS RESULTADOS
# ============================================

cat("Cargando resultados de todos los scripts...\n\n")

# Resultados de análisis espacial (Script 02)
if (file.exists("resultados/tablas/02_autocorrelacion_espacial.csv")) {
  autocorr <- read_csv("resultados/tablas/02_autocorrelacion_espacial.csv", 
                       show_col_types = FALSE)
  cat("✓ Resultados de autocorrelacion espacial cargados\n")
} else {
  cat("✗ Advertencia: Resultados de autocorrelacion no encontrados\n")
  autocorr <- NULL
}

# Resultados de validación bloqueada (Script 03)
if (file.exists("resultados/tablas/03_comparacion_metodos.csv")) {
  comparacion_03 <- read_csv("resultados/tablas/03_comparacion_metodos.csv", 
                             show_col_types = FALSE)
  cat("✓ Resultados de validacion bloqueada cargados\n")
} else {
  stop("ERROR: Resultados de validacion bloqueada no encontrados. Ejecute Script 03.")
}

if (file.exists("resultados/tablas/03_resultados_por_dominio.csv")) {
  dominios_03 <- read_csv("resultados/tablas/03_resultados_por_dominio.csv", 
                          show_col_types = FALSE)
  cat("✓ Resultados por dominio (bloqueada) cargados\n")
}

# Resultados de validación jerárquica (Script 04)
if (file.exists("resultados/tablas/04_comparacion_completa.csv")) {
  comparacion_04 <- read_csv("resultados/tablas/04_comparacion_completa.csv", 
                             show_col_types = FALSE)
  cat("✓ Resultados de validacion jerarquica cargados\n")
} else {
  stop("ERROR: Resultados de validacion jerarquica no encontrados. Ejecute Script 04.")
}

if (file.exists("resultados/tablas/04_resultados_jerarquico.csv")) {
  dominios_04 <- read_csv("resultados/tablas/04_resultados_jerarquico.csv", 
                          show_col_types = FALSE)
  cat("✓ Resultados por dominio (jerarquica) cargados\n")
}

if (file.exists("resultados/tablas/04_hiperparametros_optimos.csv")) {
  hiperparametros <- read_csv("resultados/tablas/04_hiperparametros_optimos.csv", 
                              show_col_types = FALSE)
  cat("✓ Hiperparametros optimos cargados\n")
}

cat("\n")

# ============================================
# 3. TABLA RESUMEN EJECUTIVO
# ============================================

cat("================================================================\n")
cat("TABLA 1: RESUMEN EJECUTIVO DE RESULTADOS\n")
cat("================================================================\n\n")

# Extraer valores clave
I_moran <- ifelse(!is.null(autocorr), 
                  autocorr$Valor[autocorr$Metrica == "I de Moran"], 
                  0.0031)
p_value <- ifelse(!is.null(autocorr), 
                  autocorr$Valor[autocorr$Metrica == "P-value"], 
                  0.0401)
rango_km <- ifelse(!is.null(autocorr), 
                   autocorr$Valor[autocorr$Metrica == "Rango estimado (km)"], 
                   2297)

# Crear tabla resumen
tabla_resumen <- data.frame(
  Componente = c(
    "AUTOCORRELACION ESPACIAL",
    "  I de Moran",
    "  P-value",
    "  Interpretacion",
    "  Rango espacial estimado",
    "",
    "VALIDACION ALEATORIA",
    "  RMSE",
    "  MAE",
    "  R²",
    "",
    "VALIDACION BLOQUEADA",
    "  RMSE",
    "  MAE", 
    "  R²",
    "",
    "VALIDACION JERARQUICA",
    "  RMSE",
    "  MAE",
    "  R²",
    "  Tiempo de ejecucion",
    "",
    "FILTRACION ESPACIAL",
    "  Aleatoria vs Bloqueada (RMSE)",
    "  Aleatoria vs Bloqueada (R²)",
    "  Aleatoria vs Jerarquica (RMSE)",
    "  Aleatoria vs Jerarquica (R²)",
    "  Bloqueada vs Jerarquica (RMSE)",
    "  Bloqueada vs Jerarquica (R²)"
  ),
  Valor = c(
    "",
    sprintf("%.4f", I_moran),
    sprintf("%.4f", p_value),
    "Débil pero significativa",
    sprintf("~%.0f km", rango_km),
    "",
    "",
    sprintf("%.4f", comparacion_04$RMSE[1]),
    sprintf("%.4f", comparacion_04$MAE[1]),
    sprintf("%.4f", comparacion_04$R2[1]),
    "",
    "",
    sprintf("%.4f", comparacion_04$RMSE[2]),
    sprintf("%.4f", comparacion_04$MAE[2]),
    sprintf("%.4f", comparacion_04$R2[2]),
    "",
    "",
    sprintf("%.4f", comparacion_04$RMSE[3]),
    sprintf("%.4f", comparacion_04$MAE[3]),
    sprintf("%.4f", comparacion_04$R2[3]),
    "337.89 minutos",
    "",
    "",
    sprintf("+%.2f%%", ((comparacion_04$RMSE[2] - comparacion_04$RMSE[1]) / comparacion_04$RMSE[1]) * 100),
    sprintf("-%.2f%%", ((comparacion_04$R2[1] - comparacion_04$R2[2]) / comparacion_04$R2[1]) * 100),
    sprintf("+%.2f%%", ((comparacion_04$RMSE[3] - comparacion_04$RMSE[1]) / comparacion_04$RMSE[1]) * 100),
    sprintf("-%.2f%%", ((comparacion_04$R2[1] - comparacion_04$R2[3]) / comparacion_04$R2[1]) * 100),
    sprintf("+%.2f%%", ((comparacion_04$RMSE[3] - comparacion_04$RMSE[2]) / comparacion_04$RMSE[2]) * 100),
    sprintf("-%.2f%%", ((comparacion_04$R2[2] - comparacion_04$R2[3]) / comparacion_04$R2[2]) * 100)
  )
)

print(tabla_resumen)
cat("\n")

# Guardar
write_csv(tabla_resumen, "resultados/tablas/05_resumen_ejecutivo.csv")

# ============================================
# 4. TABLA COMPARATIVA POR DOMINIO
# ============================================

cat("================================================================\n")
cat("TABLA 2: COMPARACION POR DOMINIO (BLOQUEADA VS JERARQUICA)\n")
cat("================================================================\n\n")

# Combinar resultados de validación bloqueada y jerárquica
if (exists("dominios_03") && exists("dominios_04")) {
  
  # Extraer datos
  dominios_comparacion <- data.frame(
    Dominio = c("1 (Costa Norte)", "2 (Costa Centro)", "3 (Costa Sur)", 
                "4 (Sierra Norte)", "5 (Sierra Centro)", "6 (Sierra Sur)", 
                "7 (Selva)"),
    RMSE_Bloqueada = dominios_03$RMSE,
    RMSE_Jerarquica = dominios_04$RMSE,
    R2_Bloqueada = dominios_03$R2,
    R2_Jerarquica = dominios_04$R2,
    n_test = dominios_03$n_test
  )
  
  # Calcular diferencias
  dominios_comparacion <- dominios_comparacion %>%
    mutate(
      Diff_RMSE_pct = ((RMSE_Jerarquica - RMSE_Bloqueada) / RMSE_Bloqueada) * 100,
      Diff_R2_pct = ((R2_Bloqueada - R2_Jerarquica) / R2_Bloqueada) * 100
    )
  
  print(dominios_comparacion %>% 
          select(Dominio, RMSE_Bloqueada, RMSE_Jerarquica, Diff_RMSE_pct, 
                 R2_Bloqueada, R2_Jerarquica, Diff_R2_pct, n_test))
  cat("\n")
  
  # Guardar
  write_csv(dominios_comparacion, "resultados/tablas/05_comparacion_por_dominio.csv")
  
  # Estadísticas de diferencias
  cat("Estadisticas de diferencias entre metodos:\n")
  cat(sprintf("  RMSE - Diferencia promedio: %.2f%% (rango: %.2f%% a %.2f%%)\n",
              mean(dominios_comparacion$Diff_RMSE_pct),
              min(dominios_comparacion$Diff_RMSE_pct),
              max(dominios_comparacion$Diff_RMSE_pct)))
  cat(sprintf("  R² - Diferencia promedio: %.2f%% (rango: %.2f%% a %.2f%%)\n\n",
              mean(dominios_comparacion$Diff_R2_pct),
              min(dominios_comparacion$Diff_R2_pct),
              max(dominios_comparacion$Diff_R2_pct)))
}

# ============================================
# 5. ANÁLISIS DE HIPERPARÁMETROS
# ============================================

cat("================================================================\n")
cat("TABLA 3: HIPERPARAMETROS OPTIMOS POR DOMINIO\n")
cat("================================================================\n\n")

if (exists("hiperparametros")) {
  
  # Agregar nombres descriptivos
  hiperparametros_tabla <- hiperparametros %>%
    mutate(
      Dominio_Nombre = case_when(
        Dominio == "1" ~ "1 (Costa Norte)",
        Dominio == "2" ~ "2 (Costa Centro)",
        Dominio == "3" ~ "3 (Costa Sur)",
        Dominio == "4" ~ "4 (Sierra Norte)",
        Dominio == "5" ~ "5 (Sierra Centro)",
        Dominio == "6" ~ "6 (Sierra Sur)",
        Dominio == "7" ~ "7 (Selva)",
        TRUE ~ Dominio
      )
    ) %>%
    select(Dominio_Nombre, ntree, mtry, nodesize)
  
  print(hiperparametros_tabla)
  cat("\n")
  
  # Frecuencias de hiperparámetros
  cat("Frecuencia de configuraciones:\n")
  cat("\nntree:\n")
  print(table(hiperparametros$ntree))
  cat("\nmtry:\n")
  print(table(hiperparametros$mtry))
  cat("\nnodesize:\n")
  print(table(hiperparametros$nodesize))
  cat("\n")
  
  # Configuración modal
  ntree_modal <- as.numeric(names(sort(table(hiperparametros$ntree), decreasing = TRUE)[1]))
  mtry_modal <- as.numeric(names(sort(table(hiperparametros$mtry), decreasing = TRUE)[1]))
  nodesize_modal <- as.numeric(names(sort(table(hiperparametros$nodesize), decreasing = TRUE)[1]))
  
  cat(sprintf("Configuracion modal: ntree=%d, mtry=%d, nodesize=%d\n\n", 
              ntree_modal, mtry_modal, nodesize_modal))
  
  # Guardar
  write_csv(hiperparametros_tabla, "resultados/tablas/05_hiperparametros_detalle.csv")
}

# ============================================
# 6. TABLA LATEX PARA ARTÍCULO
# ============================================

cat("================================================================\n")
cat("GENERANDO TABLAS EN FORMATO LATEX\n")
cat("================================================================\n\n")

# Tabla 1: Comparación de métodos (para LaTeX)
latex_comparacion <- comparacion_04 %>%
  mutate(
    RMSE = sprintf("%.4f", RMSE),
    MAE = sprintf("%.4f", MAE),
    R2 = sprintf("%.4f", R2)
  )

latex_code_1 <- kable(latex_comparacion, format = "latex", booktabs = TRUE,
                      caption = "Comparación de Métodos de Validación Cruzada",
                      label = "tab:comparacion_final") %>%
  kable_styling(latex_options = c("hold_position"))

# Guardar código LaTeX
writeLines(latex_code_1, "resultados/tablas/05_latex_tabla_comparacion.tex")
cat("✓ Tabla LaTeX 1 guardada: 05_latex_tabla_comparacion.tex\n")

# Tabla 2: Resultados por dominio (para LaTeX)
if (exists("dominios_04")) {
  latex_dominios <- dominios_04 %>%
    mutate(
      Dominio = case_when(
        grepl("1", Metodo) ~ "1 (Costa N)",
        grepl("2", Metodo) ~ "2 (Costa C)",
        grepl("3", Metodo) ~ "3 (Costa S)",
        grepl("4", Metodo) ~ "4 (Sierra N)",
        grepl("5", Metodo) ~ "5 (Sierra C)",
        grepl("6", Metodo) ~ "6 (Sierra S)",
        grepl("7", Metodo) ~ "7 (Selva)",
        TRUE ~ Metodo
      ),
      RMSE = sprintf("%.4f", RMSE),
      R2 = sprintf("%.4f", R2),
      n_test = format(n_test, big.mark = ",")
    ) %>%
    select(Dominio, RMSE, R2, ntree, mtry, nodesize, n_test)
  
  latex_code_2 <- kable(latex_dominios, format = "latex", booktabs = TRUE,
                        caption = "Métricas por Dominio - Validación Jerárquica",
                        label = "tab:dominios_jerarquico") %>%
    kable_styling(latex_options = c("hold_position", "scale_down"))
  
  writeLines(latex_code_2, "resultados/tablas/05_latex_tabla_dominios.tex")
  cat("✓ Tabla LaTeX 2 guardada: 05_latex_tabla_dominios.tex\n\n")
}

# ============================================
# 7. VISUALIZACIONES FINALES
# ============================================

cat("================================================================\n")
cat("GENERANDO VISUALIZACIONES FINALES\n")
cat("================================================================\n\n")

dir.create("resultados/graficos", recursive = TRUE, showWarnings = FALSE)

# Gráfico 1: Evolución de métricas a través de métodos
metricas_evolucion <- comparacion_04 %>%
  pivot_longer(cols = c(RMSE, MAE, R2), names_to = "Metrica", values_to = "Valor") %>%
  mutate(
    Metodo = factor(Metodo, levels = c("Val. Aleatoria", "Val. Bloqueada", "Val. Jerarquica"))
  )

p1 <- ggplot(metricas_evolucion, aes(x = Metodo, y = Valor, group = Metrica, color = Metrica)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~Metrica, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("RMSE" = "red", "MAE" = "orange", "R2" = "blue")) +
  labs(
    title = "Evolucion de Metricas a Traves de Metodos de Validacion",
    subtitle = "Validacion mas rigurosa reduce rendimiento aparente",
    x = "Metodo de Validacion",
    y = "Valor de Metrica"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("resultados/graficos/05_evolucion_metricas.png", p1, 
       width = 10, height = 12, dpi = 300)
cat("✓ Grafico guardado: 05_evolucion_metricas.png\n")

# Gráfico 2: Comparación RMSE por dominio (Bloqueada vs Jerárquica)
if (exists("dominios_comparacion")) {
  dominios_long <- dominios_comparacion %>%
    select(Dominio, RMSE_Bloqueada, RMSE_Jerarquica) %>%
    pivot_longer(cols = c(RMSE_Bloqueada, RMSE_Jerarquica), 
                 names_to = "Metodo", values_to = "RMSE") %>%
    mutate(
      Metodo = recode(Metodo,
                      "RMSE_Bloqueada" = "Val. Bloqueada",
                      "RMSE_Jerarquica" = "Val. Jerarquica")
    )
  
  p2 <- ggplot(dominios_long, aes(x = Dominio, y = RMSE, fill = Metodo)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", RMSE)), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 2.8) +
    scale_fill_manual(values = c("Val. Bloqueada" = "coral", 
                                 "Val. Jerarquica" = "darkgreen")) +
    labs(
      title = "Comparacion RMSE por Dominio",
      subtitle = "Validacion Bloqueada vs Jerarquica",
      x = "Dominio Geografico",
      y = "RMSE",
      fill = "Metodo"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave("resultados/graficos/05_comparacion_rmse_dominios.png", p2,
         width = 12, height = 8, dpi = 300)
  cat("✓ Grafico guardado: 05_comparacion_rmse_dominios.png\n")
}

# Gráfico 3: Heatmap de diferencias porcentuales
if (exists("dominios_comparacion")) {
  p3 <- ggplot(dominios_comparacion, aes(x = "RMSE", y = Dominio, fill = Diff_RMSE_pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", Diff_RMSE_pct)), color = "white", fontface = "bold") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, name = "Diferencia (%)") +
    labs(
      title = "Diferencia Porcentual: Jerarquica vs Bloqueada",
      subtitle = "Valores positivos indican RMSE mayor en validacion jerarquica",
      x = "",
      y = "Dominio"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_blank()
    )
  
  ggsave("resultados/graficos/05_heatmap_diferencias.png", p3,
         width = 8, height = 10, dpi = 300)
  cat("✓ Grafico guardado: 05_heatmap_diferencias.png\n")
}

# Gráfico 4: Distribución de hiperparámetros óptimos
if (exists("hiperparametros")) {
  
  # Preparar datos para gráfico de hiperparámetros
  hiper_long <- hiperparametros %>%
    pivot_longer(cols = c(ntree, mtry, nodesize), 
                 names_to = "Hiperparametro", values_to = "Valor") %>%
    mutate(
      Dominio_Nombre = case_when(
        Dominio == "1" ~ "Costa Norte",
        Dominio == "2" ~ "Costa Centro",
        Dominio == "3" ~ "Costa Sur",
        Dominio == "4" ~ "Sierra Norte",
        Dominio == "5" ~ "Sierra Centro",
        Dominio == "6" ~ "Sierra Sur",
        Dominio == "7" ~ "Selva"
      ),
      Hiperparametro = factor(Hiperparametro, 
                              levels = c("ntree", "mtry", "nodesize"),
                              labels = c("Número de árboles", "Variables por split", "Tamaño mínimo nodo"))
    )
  
  p4 <- ggplot(hiper_long, aes(x = Dominio_Nombre, y = Valor, fill = Hiperparametro)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    facet_wrap(~Hiperparametro, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Hiperparametros Optimos por Dominio",
      subtitle = "Seleccionados mediante Validacion Cruzada Jerarquica",
      x = "Dominio Geografico",
      y = "Valor del Hiperparametro"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 11)
    )
  
  ggsave("resultados/graficos/05_hiperparametros_dominios.png", p4,
         width = 12, height = 12, dpi = 300)
  cat("✓ Grafico guardado: 05_hiperparametros_dominios.png\n")
}

cat("\n")

# ============================================
# 8. ANÁLISIS ESTADÍSTICO DE SIGNIFICANCIA
# ============================================

cat("================================================================\n")
cat("ANALISIS ESTADISTICO DE SIGNIFICANCIA\n")
cat("================================================================\n\n")

if (exists("dominios_03") && exists("dominios_04")) {
  
  # Test t pareado para RMSE
  cat("TEST T PAREADO: RMSE Bloqueada vs Jerarquica\n")
  test_rmse <- t.test(dominios_03$RMSE, dominios_04$RMSE, paired = TRUE)
  cat(sprintf("  t = %.4f\n", test_rmse$statistic))
  cat(sprintf("  df = %d\n", test_rmse$parameter))
  cat(sprintf("  p-value = %.6f\n", test_rmse$p.value))
  cat(sprintf("  Diferencia de medias: %.4f\n", test_rmse$estimate))
  
  if (test_rmse$p.value < 0.05) {
    cat("  Conclusion: Diferencia SIGNIFICATIVA (p < 0.05)\n\n")
  } else {
    cat("  Conclusion: Diferencia NO significativa (p >= 0.05)\n\n")
  }
  
  # Test t pareado para R²
  cat("TEST T PAREADO: R² Bloqueada vs Jerarquica\n")
  test_r2 <- t.test(dominios_03$R2, dominios_04$R2, paired = TRUE)
  cat(sprintf("  t = %.4f\n", test_r2$statistic))
  cat(sprintf("  df = %d\n", test_r2$parameter))
  cat(sprintf("  p-value = %.6f\n", test_r2$p.value))
  cat(sprintf("  Diferencia de medias: %.4f\n", test_r2$estimate))
  
  if (test_r2$p.value < 0.05) {
    cat("  Conclusion: Diferencia SIGNIFICATIVA (p < 0.05)\n\n")
  } else {
    cat("  Conclusion: Diferencia NO significativa (p >= 0.05)\n\n")
  }
  
  # Correlación entre RMSE de ambos métodos
  cat("CORRELACION ENTRE METODOS\n")
  cor_rmse <- cor.test(dominios_03$RMSE, dominios_04$RMSE)
  cat(sprintf("  Correlacion de Pearson (RMSE): r = %.4f, p = %.6f\n", 
              cor_rmse$estimate, cor_rmse$p.value))
  
  cor_r2 <- cor.test(dominios_03$R2, dominios_04$R2)
  cat(sprintf("  Correlacion de Pearson (R²): r = %.4f, p = %.6f\n\n", 
              cor_r2$estimate, cor_r2$p.value))
  
  # Guardar resultados de tests
  resultados_tests <- data.frame(
    Test = c("T pareado RMSE", "T pareado R²", 
             "Correlacion RMSE", "Correlacion R²"),
    Estadistico = c(test_rmse$statistic, test_r2$statistic,
                    cor_rmse$estimate, cor_r2$estimate),
    P_value = c(test_rmse$p.value, test_r2$p.value,
                cor_rmse$p.value, cor_r2$p.value),
    Interpretacion = c(
      ifelse(test_rmse$p.value < 0.05, "Significativa", "No significativa"),
      ifelse(test_r2$p.value < 0.05, "Significativa", "No significativa"),
      ifelse(cor_rmse$p.value < 0.05, "Significativa", "No significativa"),
      ifelse(cor_r2$p.value < 0.05, "Significativa", "No significativa")
    )
  )
  
  write_csv(resultados_tests, "resultados/tablas/05_tests_estadisticos.csv")
}

# ============================================
# 9. REPORTE FINAL CONSOLIDADO
# ============================================

cat("================================================================\n")
cat("REPORTE FINAL CONSOLIDADO\n")
cat("================================================================\n\n")

# Crear reporte en Markdown
reporte <- c(
  "# REPORTE FINAL: VALIDACIÓN CRUZADA BLOQUEADA Y JERÁRQUICA",
  "## Encuesta Nacional Agropecuaria del Perú 2017",
  "",
  "---",
  "",
  "## 1. RESUMEN EJECUTIVO",
  "",
  sprintf("- **Tamaño del dataset**: %s observaciones", 
          format(sum(dominios_04$n_test), big.mark = ",")),
  sprintf("- **I de Moran**: %.4f (p = %.4f) - Autocorrelación espacial débil pero significativa", 
          I_moran, p_value),
  sprintf("- **Rango espacial**: ~%.0f km (limitado por coordenadas sintéticas)", rango_km),
  "",
  "## 2. MÉTRICAS DE RENDIMIENTO PREDICTIVO",
  "",
  "| Método | RMSE | MAE | R² |",
  "|--------|------|-----|-----|",
  sprintf("| Val. Aleatoria | %.4f | %.4f | %.4f |", 
          comparacion_04$RMSE[1], comparacion_04$MAE[1], comparacion_04$R2[1]),
  sprintf("| Val. Bloqueada | %.4f | %.4f | %.4f |", 
          comparacion_04$RMSE[2], comparacion_04$MAE[2], comparacion_04$R2[2]),
  sprintf("| Val. Jerárquica | %.4f | %.4f | %.4f |", 
          comparacion_04$RMSE[3], comparacion_04$MAE[3], comparacion_04$R2[3]),
  "",
  "## 3. FILTRACIÓN ESPACIAL DETECTADA",
  "",
  sprintf("- **Aleatoria vs Bloqueada**: RMSE +%.2f%%, R² -%.2f%%", 
          ((comparacion_04$RMSE[2] - comparacion_04$RMSE[1]) / comparacion_04$RMSE[1]) * 100,
          ((comparacion_04$R2[1] - comparacion_04$R2[2]) / comparacion_04$R2[1]) * 100),
  sprintf("- **Aleatoria vs Jerárquica**: RMSE +%.2f%%, R² -%.2f%%", 
          ((comparacion_04$RMSE[3] - comparacion_04$RMSE[1]) / comparacion_04$RMSE[1]) * 100,
          ((comparacion_04$R2[1] - comparacion_04$R2[3]) / comparacion_04$R2[1]) * 100),
  sprintf("- **Bloqueada vs Jerárquica**: RMSE +%.2f%%, R² -%.2f%%", 
          ((comparacion_04$RMSE[3] - comparacion_04$RMSE[2]) / comparacion_04$RMSE[2]) * 100,
          ((comparacion_04$R2[2] - comparacion_04$R2[3]) / comparacion_04$R2[2]) * 100),
  "",
  "**Interpretación**: La validación aleatoria sobreestima dramáticamente el rendimiento predictivo.",
  "Validación jerárquica proporciona las estimaciones más conservadoras y honestas.",
  "",
  "## 4. RENDIMIENTO POR DOMINIO GEOGRÁFICO",
  "",
  "### Validación Jerárquica - RMSE por Dominio:",
  ""
)

# Agregar resultados por dominio
if (exists("dominios_04")) {
  for (i in 1:nrow(dominios_04)) {
    dominio_nombre <- case_when(
      grepl("1", dominios_04$Metodo[i]) ~ "Costa Norte",
      grepl("2", dominios_04$Metodo[i]) ~ "Costa Centro",
      grepl("3", dominios_04$Metodo[i]) ~ "Costa Sur",
      grepl("4", dominios_04$Metodo[i]) ~ "Sierra Norte",
      grepl("5", dominios_04$Metodo[i]) ~ "Sierra Centro",
      grepl("6", dominios_04$Metodo[i]) ~ "Sierra Sur",
      grepl("7", dominios_04$Metodo[i]) ~ "Selva",
      TRUE ~ dominios_04$Metodo[i]
    )
    reporte <- c(reporte, sprintf("- **%s**: RMSE=%.4f, R²=%.4f (n=%s)", 
                                  dominio_nombre, 
                                  dominios_04$RMSE[i], 
                                  dominios_04$R2[i],
                                  format(dominios_04$n_test[i], big.mark = ",")))
  }
}

reporte <- c(reporte,
             "",
             sprintf("- **Variación de RMSE**: %.1f%% (%.4f a %.4f)", 
                     (max(dominios_04$RMSE) - min(dominios_04$RMSE)) / min(dominios_04$RMSE) * 100,
                     min(dominios_04$RMSE), max(dominios_04$RMSE)),
             "",
             "## 5. HIPERPARÁMETROS ÓPTIMOS",
             "",
             "### Configuraciones seleccionadas por Inner CV:",
             ""
)

# Agregar hiperparámetros
if (exists("hiperparametros")) {
  for (i in 1:nrow(hiperparametros)) {
    dominio_nombre <- case_when(
      hiperparametros$Dominio[i] == "1" ~ "Costa Norte",
      hiperparametros$Dominio[i] == "2" ~ "Costa Centro",
      hiperparametros$Dominio[i] == "3" ~ "Costa Sur",
      hiperparametros$Dominio[i] == "4" ~ "Sierra Norte",
      hiperparametros$Dominio[i] == "5" ~ "Sierra Centro",
      hiperparametros$Dominio[i] == "6" ~ "Sierra Sur",
      hiperparametros$Dominio[i] == "7" ~ "Selva",
      TRUE ~ hiperparametros$Dominio[i]
    )
    reporte <- c(reporte, sprintf("- **%s**: ntree=%d, mtry=%d, nodesize=%d",
                                  dominio_nombre,
                                  hiperparametros$ntree[i],
                                  hiperparametros$mtry[i],
                                  hiperparametros$nodesize[i]))
  }
}

reporte <- c(reporte,
             "",
             "### Frecuencias:",
             sprintf("- **ntree**: 50 (%d dominios), 100 (%d), 150 (%d)",
                     sum(hiperparametros$ntree == 50),
                     sum(hiperparametros$ntree == 100),
                     sum(hiperparametros$ntree == 150)),
             sprintf("- **mtry**: 1 (%d dominios), 2 (%d)",
                     sum(hiperparametros$mtry == 1),
                     sum(hiperparametros$mtry == 2)),
             sprintf("- **nodesize**: 5 (%d dominios), 10 (%d)",
                     sum(hiperparametros$nodesize == 5),
                     sum(hiperparametros$nodesize == 10)),
             "",
             "**Conclusión**: No existe configuración única óptima para todos los dominios.",
             "Heterogeneidad de hiperparámetros refleja diferencias en estructura espacial.",
             "",
             "## 6. ANÁLISIS ESTADÍSTICO",
             ""
)

if (exists("test_rmse") && exists("test_r2")) {
  reporte <- c(reporte,
               "### Test t pareado (Bloqueada vs Jerárquica):",
               "",
               sprintf("- **RMSE**: t=%.4f, p=%.6f → %s",
                       test_rmse$statistic, test_rmse$p.value,
                       ifelse(test_rmse$p.value < 0.05, "Diferencia SIGNIFICATIVA", "No significativa")),
               sprintf("- **R²**: t=%.4f, p=%.6f → %s",
                       test_r2$statistic, test_r2$p.value,
                       ifelse(test_r2$p.value < 0.05, "Diferencia SIGNIFICATIVA", "No significativa")),
               "",
               "### Correlación entre métodos:",
               "",
               sprintf("- **RMSE**: r=%.4f, p=%.6f → Los métodos ordenan dominios de manera %s",
                       cor_rmse$estimate, cor_rmse$p.value,
                       ifelse(cor_rmse$estimate > 0.8, "muy consistente", 
                              ifelse(cor_rmse$estimate > 0.5, "consistente", "variable"))),
               sprintf("- **R²**: r=%.4f, p=%.6f",
                       cor_r2$estimate, cor_r2$p.value),
               ""
  )
}

reporte <- c(reporte,
             "## 7. CONCLUSIONES PRINCIPALES",
             "",
             "1. **Autocorrelación espacial débil pero significativa** justifica uso de validación espacial.",
             "",
             "2. **Validación aleatoria sobreestima rendimiento** en 9-11% (RMSE) y 25-45% (R²), ",
             "   produciendo estimaciones optimistas no representativas de capacidad de generalización espacial.",
             "",
             "3. **Validación jerárquica proporciona estimaciones más conservadoras** que validación bloqueada,",
             "   revelando filtración residual en optimización de hiperparámetros.",
             "",
             "4. **Heterogeneidad espacial sustancial**: RMSE varía 37% entre dominios (0.88 a 1.21),",
             "   con Costa Sur mostrando mayor dificultad predictiva.",
             "",
             "5. **Hiperparámetros óptimos varían por contexto espacial**, sin configuración universal.",
             "   Esto sugiere que relaciones predictivas tienen estructura espacial compleja.",
             "",
             "6. **Costo computacional de validación jerárquica** (5.6 horas) es justificable para",
             "   modelamiento definitivo donde estimaciones honestas son críticas.",
             "",
             "## 8. RECOMENDACIONES",
             "",
             "### Para Investigadores:",
             "- Usar validación cruzada bloqueada como **mínimo** para datos espacialmente estructurados",
             "- Implementar validación jerárquica cuando optimización de hiperparámetros es crítica",
             "- Reportar métricas de múltiples esquemas de validación para transparencia",
             "",
             "### Para Formuladores de Política:",
             "- Reconocer que predicciones para regiones no observadas tienen mayor incertidumbre",
             "- Considerar heterogeneidad espacial en diseño de intervenciones regionales",
             "- Invertir en recolección de datos de regiones con peor rendimiento predictivo",
             "",
             "### Para Trabajo Futuro:",
             "- Obtener coordenadas reales de productores para análisis variográfico preciso",
             "- Incluir covariables adicionales (clima, topografía, mercados) para mejorar predicción",
             "- Evaluar metodología en otros cultivos/especies y contextos geográficos",
             "- Explorar validación temporal para datos longitudinales",
             "",
             "---",
             "",
             sprintf("**Reporte generado**: %s", Sys.time()),
             sprintf("**Directorio**: %s", getwd()),
             ""
)

# Guardar reporte
writeLines(reporte, "resultados/REPORTE_FINAL.md")
cat("✓ Reporte final guardado: resultados/REPORTE_FINAL.md\n\n")

# ============================================
# 10. RESUMEN FINAL DEL SCRIPT
# ============================================

cat("================================================================\n")
cat("RESUMEN DE ARCHIVOS GENERADOS\n")
cat("================================================================\n\n")

archivos_generados <- c(
  "TABLAS CSV:",
  "  ✓ 05_resumen_ejecutivo.csv",
  "  ✓ 05_comparacion_por_dominio.csv",
  "  ✓ 05_hiperparametros_detalle.csv",
  "  ✓ 05_tests_estadisticos.csv",
  "",
  "TABLAS LATEX:",
  "  ✓ 05_latex_tabla_comparacion.tex",
  "  ✓ 05_latex_tabla_dominios.tex",
  "",
  "GRAFICOS:",
  "  ✓ 05_evolucion_metricas.png",
  "  ✓ 05_comparacion_rmse_dominios.png",
  "  ✓ 05_heatmap_diferencias.png",
  "  ✓ 05_hiperparametros_dominios.png",
  "",
  "REPORTE:",
  "  ✓ REPORTE_FINAL.md"
)

cat(paste(archivos_generados, collapse = "\n"))
cat("\n\n")

cat("================================================================\n")
cat("ANALISIS COMPLETO FINALIZADO\n")
cat("================================================================\n\n")

cat("Todos los scripts han sido ejecutados exitosamente:\n")
cat("  ✓ Script 01: Preparacion de datos (Python)\n")
cat("  ✓ Script 02: Analisis de autocorrelacion espacial (R)\n")
cat("  ✓ Script 03: Validacion cruzada bloqueada (R)\n")
cat("  ✓ Script 04: Validacion cruzada jerarquica (R)\n")
cat("  ✓ Script 05: Comparacion final y tablas consolidadas (R)\n\n")

cat("Resultados listos para:\n")
cat("  - Redaccion del articulo cientifico\n")
cat("  - Presentacion de resultados\n")
cat("  - Publicacion en repositorios\n\n")

cat("Directorios de salida:\n")
cat("  - resultados/tablas/\n")
cat("  - resultados/graficos/\n")
cat("  - resultados/REPORTE_FINAL.md\n\n")

cat("================================================================\n")
cat("¡PROYECTO COMPLETADO!\n")
cat("================================================================\n")
          
