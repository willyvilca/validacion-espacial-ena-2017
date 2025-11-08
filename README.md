# Validación Cruzada Espacial para Datos Agropecuarios del Perú

[![DOI](https://img.shields.io/badge/DOI-pendiente-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 📄 Resumen

Código reproducible del artículo "Validación Cruzada Bloqueada y Jerárquica para Datos Agropecuarios con Autocorrelación Espacial: Aplicación a la Encuesta Nacional Agropecuaria del Perú".

**Autores:** Willy Vilca Apaza, Fred Torres Cruz, Edgar Eloy Carpio Vargas  
**Institución:** Universidad Nacional del Altiplano, Puno, Perú

## 🎯 Objetivos

1. Cuantificar autocorrelación espacial en inventario ganadero
2. Comparar validación cruzada aleatoria vs. bloqueada
3. Evaluar costo-beneficio de validación jerárquica

## 📊 Resultados Principales

- Autocorrelación espacial débil pero significativa (I=0.0031, p=0.0401)
- Validación aleatoria sobreestima rendimiento en 9.1% (RMSE) y 25.1% (R²)
- Validación jerárquica 60× más costosa sin beneficio estadístico significativo

## 🗂️ Estructura del Repositorio
```
├── articulo/      # Artículo en LaTeX y PDF
├── code/          # Scripts R reproducibles
├── data/          # Datos procesados (ENA 2017)
├── figuras/       # Figuras del artículo
└── resultados/    # Tablas y modelos guardados
```

## 🚀 Reproducibilidad

### Requisitos

- R >= 4.3.0
- RStudio (recomendado)

### Paquetes R necesarios
```r
install.packages(c("randomForest", "spdep", "gstat", "caret", 
                   "ggplot2", "dplyr", "sf"))
```

### Ejecución
```r
# 1. Limpiar datos
source("code/01_limpieza_datos.R")

# 2. Análisis espacial
source("code/02_analisis_espacial.R")

# 3. Validación cruzada
source("code/03_validacion_cruzada.R")

# 4. Generar figuras
source("code/04_visualizaciones.R")
```

## 📥 Datos

Los microdatos de ENA 2017 son de acceso público mediante solicitud a:
- **INEI:** https://www.inei.gob.pe

⚠️ Por confidencialidad, no se incluyen coordenadas geográficas reales.

## 📖 Citación
```bibtex
@inproceedings{vilca2025validacion,
  title={Validación Cruzada Bloqueada y Jerárquica para Datos Agropecuarios},
  author={Vilca Apaza, Willy and Torres Cruz, Fred and Carpio Vargas, Edgar E.},
  booktitle={Proceedings of [Nombre Conferencia]},
  year={2025},
  publisher={Springer}
}
```

## 📧 Contacto

- Willy Vilca Apaza: w.vilca@unap.edu.pe
- Universidad Nacional del Altiplano, Puno, Perú

## 📜 Licencia

Este proyecto está bajo licencia MIT. Ver [LICENSE](LICENSE) para detalles.