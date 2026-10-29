# Proyecto_DC
Proyecto final data science
# ETc del Cultivo por FAO-56 — Doble Kc y Riego al 50 % TAW

Proyecto de evapotranspiración del cultivo (ETc) desarrollado en Python, aplicando el método **FAO-56** (Allen et al., 1998) para un único escenario con riego.

## 📋 Descripción

Este proyecto calcula la evapotranspiración del cultivo usando:
- **ETo** por Penman-Monteith FAO-56 (Eq. 6).
- **Coeficiente doble**: $Kc = Kcb + Ke$ (basal + evaporación del suelo).
- **ETc en condiciones no estándar** con coeficiente de estrés hídrico $Ks$.
- **Riego**: reposición a capacidad de campo (FC) cuando la depleción de la zona radicular alcanza el **50 % de TAW**.

## 🌱 Datos de entrada

Se utilizaron datos climáticos de la hoja **"ETc no estándar"** de `Datos_clima_P2.xlsx`:

| Variable | Descripción |
|----------|-------------|
| Periodo | 15-09-2025 al 07-04-2026 (205 días) |
| Tmin / Tmax | Temperatura mínima y máxima del aire [°C] |
| RHmin / RHmax | Humedad relativa mínima y máxima [%] |
| P | Precipitación [mm] |
| Rs | Radiación solar [MJ m⁻² día⁻¹] |
| u2 | Velocidad del viento a 2 m [km/h] |
| ETo | Evapotranspiración de referencia [mm/día] |

Sitio: latitud -33.45°, altitud 520 msnm.

## 🧮 Metodología

### 1. ETo (Penman-Monteith FAO-56)
Se calcula la evaporación de referencia diaria a partir de temperatura, humedad, radiación solar y viento.

### 2. Coeficiente doble Kc
- **Kcb**: coeficiente basal de cultivo, varía por etapas fenológicas (inicio, desarrollo, medio, final).
- **Ke**: coeficiente de evaporación del suelo, depende de la depleción de la capa superficial ($De$).
- **Kc = Kcb + Ke**

### 3. ETc en condiciones no estándar
- **ETc estándar** = $(Kcb + Ke) \cdot ETo$ (sin estrés hídrico)
- **ETc no estándar** = $(Ks \cdot Kcb + Ke) \cdot ETo$ (con estrés hídrico)

Donde $Ks$ es el coeficiente de estrés hídrico (Eq. 82 FAO-56), calculado a partir de la depleción de la zona radicular ($Dr$).

### 4. Condición de riego
- Umbral: **50 % de TAW** (agua total disponible en la zona radicular).
- Al alcanzar el umbral, se aplica riego de reposición a capacidad de campo ($Dr \to 0$, $De \to 0$).

## 📊 Resultados principales (escenario con riego)

| Indicador | Valor |
|-----------|-------|
| ETo media | 3.87 mm/día |
| ETc estándar total | 791.3 mm |
| ETc no estándar total | 791.3 mm |
| Riego total aplicado | 684.3 mm |
| Número de riegos | 11 eventos |
| Días con estrés (Ks < 0.95) | 0 de 205 |

Como el riego repone el suelo al 50 % de TAW, el cultivo no sufre estrés hídrico ($Ks \approx 1$ todo el ciclo), por lo que ETc estándar ≈ ETc no estándar.

## 📁 Estructura del repositorio
