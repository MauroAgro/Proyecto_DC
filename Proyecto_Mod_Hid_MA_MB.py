"""
Proyecto_Mod_Hid_MA_MB.py

=========================

Cálculo de la evapotranspiración del cultivo (ETc) por el método FAO-56

(Allen et al., 1998, "Crop evapotranspiration - FAO Irrigation and Drainage
Paper 56"), para UN SOLO escenario, considerando:

  1. ETo por la ecuación de Penman-Monteith FAO-56.
  2. Coeficiente DUAL de cultivo: Kc = Kcb + Ke
        Kcb = coeficiente basal de cultivo (transpiración)
        Ke  = coeficiente de evaporación del suelo
  3. ETc en condiciones NO ESTÁNDAR con coeficiente de estrés hídrico Ks:
        ETc_estándar    = (Kcb + Ke) · ETo          (sin limitación hídrica)
        ETc_no_estándar = (Ks · Kcb + Ke) · ETo     (con estrés hídrico)

CONDICIÓN DE RIEGO (este escenario)
-----------------------------------
  - Riego por reposición a capacidad de campo (FC).
  - Umbral de riego: se riega cuando la depleción de la zona radicular (Dr)
    alcanza el 50 % del TAW (agua total disponible).
  - Al regar, se restaura Dr a 0 (zona radicular a FC) y se asume que el riego
    también humedece la capa superficial (De → 0).

Se incluye un balance hídrico diario simplificado del suelo que actualiza:
  - De: depleción de agua de la capa superficial (controla Ke)  [Eq. 77-79 FAO-56]
  - Dr: depleción de agua en la zona radicular (controla Ks)    [Eq. 81-84 FAO-56]

Datos de entrada
----------------
Proyecto_Final/Datos_clima_P2.xlsx, hoja "ETc no estándar":
  periodo 15-09-2025 al 07-04-2026 (205 días)
  columnas: fecha, Tmin, Tmax, RHmin, RHmax, precipitación, radiación solar,
            velocidad de viento, ETo (de referencia, solo para comparación).

Requisitos
----------
numpy, pandas, matplotlib

Autor: proyecto de evapotranspiración del cultivo (FAO-56)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ===========================================================================
# PARÁMETROS DEL SITIO (ajustables)
# ===========================================================================
LAT = -33.45       # latitud [grados decimales; negativo = hemisferio sur]
Z   = 520.0        # altitud [msnm]
ALBEDO = 0.23      # albedo de la superficie de referencia (césped FAO-56)

# ===========================================================================
# CULTIVO Y SUELO (representativos; ajustar según el cultivo real)
# ===========================================================================
@dataclass
class CropSoil:
    # --- Fenología (días). Suma = 205, coincide con el periodo de datos ---
    L_ini: int = 25   # etapa inicial
    L_dev: int = 45   # desarrollo
    L_mid: int = 90   # estación media
    L_end: int = 45   # final de estación

    # --- Coeficiente basal de cultivo Kcb (doble Kc, aprox. Tabla 17 FAO-56) ---
    Kcb_ini: float = 0.15
    Kcb_mid: float = 1.05
    Kcb_end: float = 0.90

    # --- Evaporación máxima del suelo con superficie mojada ---
    Ke_max: float = 0.12   # valor representativo (FAO-56 da un rango según clima)

    # --- Suelo ---
    Zr: float = 0.80       # profundidad de la zona radicular [m]
    th_FC: float = 0.28    # humedad en capacidad de campo [m3/m3]
    th_WP: float = 0.13    # punto de marchitez permanente [m3/m3]
    th_01: float = 0.10    # humedad a -0.1 bar (para TEW) [m3/m3]
    Ze: float = 0.15       # profundidad de la capa evaporable del suelo [m]
    p_tab: float = 0.30    # fracción de TAW antes de que aparezca estrés (Tabla 22)

    # --- Condiciones iniciales (día 0) ---
    Dr0: float = 0.0   # depleción inicial de la zona radicular [mm] (0 = FC)
    De0: float = 0.0   # depleción inicial de la capa superficial [mm]


# ===========================================================================
# FUNCIONES FÍSICAS FAO-56
# ===========================================================================
def es_sat(T):
    """Presión de vapor de saturación es(T) [kPa] (Eq. 11 FAO-56)."""
    return 0.6108 * np.exp(17.27 * T / (T + 237.3))


def slope_vp(T):
    """Pendiente de la curva de presiones de vapor Δ [kPa/°C] (Eq. 13)."""
    es = es_sat(T)
    return 4098 * es / (T + 237.3) ** 2


def atmospheric_pressure(z):
    """Presión atmosférica P [kPa] a altitud z [m] (Eq. 7)."""
    return 101.3 * ((293.0 - 0.0065 * z) / 293.0) ** 5.26


def psychrometric_constant(P):
    """Constante psicrométrica γ [kPa/°C] (Eq. 8)."""
    return 0.665e-3 * P


def extraterrestrial_radiation(doy, lat_deg):
    """
    Radiación extraterrestre Ra [MJ m-2 d-1] (Eq. 21-31 FAO-56).
    doy: día del año (1-365); lat_deg: latitud [grados].
    """
    lat = np.radians(lat_deg)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365.0)          # distancia relativa Tierra-Sol
    delta = 0.409 * np.sin(2 * np.pi * doy / 365.0 - 1.39)    # declinación solar [rad]
    ws = np.arccos(-np.tan(lat) * np.tan(delta))              # ángulo horario del ocaso [rad]
    Gsc = 0.0820                                              # constante solar [MJ m-2 min-1]
    Ra = (24 * 60 / np.pi) * Gsc * dr * (
        ws * np.sin(lat) * np.sin(delta)
        + np.cos(lat) * np.cos(delta) * np.sin(ws)
    )
    return float(Ra)


def compute_eto(tmin, tmax, rhmin, rhmax, rs, u2_kmh, doy, lat=LAT, z=Z):
    """
    ETo por Penman-Monteith FAO-56 [mm/día] (Eq. 6).
    u2_kmh: velocidad del viento [km/h] (asumida a 2 m de altura).
    """
    T = (tmin + tmax) / 2.0
    u2 = u2_kmh / 3.6                                          # [m/s]

    Delta = slope_vp(T)                                        # [kPa/°C]
    P = atmospheric_pressure(z)                                # [kPa]
    gamma = psychrometric_constant(P)                          # [kPa/°C]

    # Presión de vapor de saturación media y real [kPa] (Eq. 11, 15, 17)
    es = (es_sat(tmax) + es_sat(tmin)) / 2.0
    ea = (es_sat(tmin) * rhmax / 100.0 +
          es_sat(tmax) * rhmin / 100.0) / 2.0
    ea = float(np.clip(ea, 0.0, es))

    # Radiación neta (G ≈ 0 para paso diario)
    Rns = (1.0 - ALBEDO) * rs                                  # [MJ m-2 d-1]
    Ra = extraterrestrial_radiation(doy, lat)                  # [MJ m-2 d-1]
    Rso = (0.75 + 2e-5 * z) * Ra                               # [Eq. 38]
    Rnl = (4.903e-9 * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2.0
           * (0.34 - 0.14 * np.sqrt(ea))
           * (1.35 * rs / Rso - 0.35))                        # [Eq. 55]
    Rnl = float(np.clip(Rnl, 0.0, None))
    Rn = Rns - Rnl                                             # [MJ m-2 d-1]

    numer = (0.408 * Delta * (Rn - 0.0)
             + gamma * (900.0 / (T + 273.0)) * u2 * (es - ea))
    denom = Delta + gamma * (1.0 + 0.34 * u2)
    eto = numer / denom
    return float(np.clip(eto, 0.0, None))


# ===========================================================================
# COEFICIENTE BASAL Kcb POR ETAPA (doble Kc)
# ===========================================================================
def kcb_of_day(i, cs):
    """
    Coeficiente basal Kcb para el día i (contado desde la siembra, 0-based),
    variando linealmente en desarrollo y final de estación.
    """
    L1 = cs.L_ini
    L2 = L1 + cs.L_dev
    L3 = L2 + cs.L_mid
    if i < L1:
        return cs.Kcb_ini
    elif i < L2:                       # desarrollo: línea Kcb_ini -> Kcb_mid
        frac = (i - L1) / cs.L_dev
        return cs.Kcb_ini + frac * (cs.Kcb_mid - cs.Kcb_ini)
    elif i < L3:                       # estación media
        return cs.Kcb_mid
    else:                              # final: línea Kcb_mid -> Kcb_end
        frac = (i - L3) / cs.L_end
        return cs.Kcb_mid + frac * (cs.Kcb_end - cs.Kcb_mid)


# ===========================================================================
# ESCENARIO ÚNICO: simulación diaria FAO-56 (doble Kc + riego al 50% TAW)
# ===========================================================================
def run_scenario(df, cs, lat=LAT, z=Z,
                 riego_umbral_frac=0.5, aplicar_riego=True):
    """
    Ejecuta la simulación día a día.

    Parámetros
    ----------
    riego_umbral_frac : float
        Fracción de TAW que dispara el riego (0.5 = 50 % de TAW).
    aplicar_riego : bool
        Si True, se aplica riego de reposición a FC al alcanzar el umbral.

    Retorna un DataFrame con:
      fecha, doy, Tmin, Tmax, RHmin, RHmax, P, Rs, u2, ETo,
      riego_mm, Kcb, Ke, Kc, Ks,
      ETc_estandar, ETc_no_estandar,
      De (depleción capa superficial), Dr (depleción zona radicular),
      ETc_estandar_cum, ETc_no_estandar_cum, riego_cum (cumuladas [mm])
    """
    n = len(df)

    # Parámetros derivados del suelo [mm]
    TAW = 1000.0 * (cs.th_FC - cs.th_WP) * cs.Zr   # agua total disponible zona raíz
    TEW = 1000.0 * (cs.th_FC - cs.th_01) * cs.Ze   # agua evaporable total capa sup.
    REW = 0.5 * TEW                                # agua evaporable fácilmente disponible
    umbral_riego = riego_umbral_frac * TAW         # umbral que dispara el riego

    # Arreglos de salida
    ETo = np.zeros(n)
    Kcb = np.zeros(n)
    Ke = np.zeros(n)
    Kc = np.zeros(n)
    Ks = np.ones(n)
    ETc_std = np.zeros(n)
    ETc_ns = np.zeros(n)
    De = np.zeros(n)
    Dr = np.zeros(n)
    riego = np.zeros(n)

    De_prev = cs.De0
    Dr_prev = cs.Dr0

    for i in range(n):
        row = df.iloc[i]
        tmin = float(row["Tmin"])
        tmax = float(row["Tmax"])
        rhmin = float(row["RHmin"])
        rhmax = float(row["RHmax"])
        prec = float(row["P"])          # precipitación [mm]
        rs = float(row["Rs"])           # radiación solar [MJ/m2/día]
        u2 = float(row["u2"])           # viento [km/h]
        doy = int(row["doy"])

        # ------------------------------------------------------------------
        # RIEGO al inicio del día: reposición a capacidad de campo
        # ------------------------------------------------------------------
        riego_i = 0.0
        if aplicar_riego and Dr_prev >= umbral_riego:
            riego_i = Dr_prev           # reponer lo faltante para llegar a FC
            Dr_start = 0.0              # zona radicular a FC (Dr = 0)
            De_start = 0.0              # capa superficial mojada por el riego
        else:
            Dr_start = Dr_prev
            De_start = De_prev

        # ------------------------------------------------------------------
        # 1. ETo (Penman-Monteith FAO-56)
        # ------------------------------------------------------------------
        eto = compute_eto(tmin, tmax, rhmin, rhmax, rs, u2, doy, lat, z)
        ETo[i] = eto

        # ------------------------------------------------------------------
        # 2. Coeficiente basal Kcb y evaporación del suelo Ke
        # ------------------------------------------------------------------
        kcb = kcb_of_day(i, cs)
        Kcb[i] = kcb

        # Ke según depleción de la capa superficial al inicio del día (Eq. 77-79)
        if De_start <= REW:
            ke = cs.Ke_max                       # superficie húmeda
        elif De_start < TEW:
            ke = cs.Ke_max * (TEW - De_start) / (TEW - REW)
        else:
            ke = 0.0                             # capa superficial seca
        Ke[i] = ke

        kc = kcb + ke
        Kc[i] = kc

        # ------------------------------------------------------------------
        # 3. Coeficiente de estrés hídrico Ks (Eq. 82, sobre Dr al inicio del día)
        # ------------------------------------------------------------------
        p = cs.p_tab + 0.04 * (eto - 5.0)        # ajuste climático Eq. 84
        p = float(np.clip(p, 0.10, 0.80))
        DRc = (1.0 - p) * TAW                    # depleción crítica
        if Dr_start <= DRc:
            ks = 1.0
        else:
            ks = (TAW - Dr_start) / (p * TAW)
            ks = float(np.clip(ks, 0.0, 1.0))
        Ks[i] = ks

        # ------------------------------------------------------------------
        # 4. ETc estándar y no estándar
        # ------------------------------------------------------------------
        etc_std = (kcb + ke) * eto               # Ks = 1 (condiciones estándar)
        etc_ns = (ks * kcb + ke) * eto           # con estrés hídrico (no estándar)
        ETc_std[i] = etc_std
        ETc_ns[i] = etc_ns

        # ------------------------------------------------------------------
        # 5. Balance hídrico diario del suelo (con precipitación efectiva)
        # ------------------------------------------------------------------
        # Capa superficial: gana agua con precipitación, pierde con Ke·ETo
        espacio_surf = TEW - De_start
        w_surf = min(prec, espacio_surf)         # agua que rellena la capa superficial
        De_new = De_start - w_surf + ke * eto
        De_new = float(np.clip(De_new, 0.0, TEW))
        De[i] = De_new

        # Zona radicular: pierde agua por transpiración real (Ks·Kcb·ETo)
        # y gana el exceso de precipitación que supera la capa superficial.
        w_to_root = max(0.0, prec - w_surf)
        Dr_new = Dr_start - w_to_root + (ks * kcb * eto)
        Dr_new = float(np.clip(Dr_new, 0.0, TAW))   # exceso = drenaje profundo (perdido)
        Dr[i] = Dr_new

        riego[i] = riego_i

        # Preparar siguiente día
        De_prev = De_new
        Dr_prev = Dr_new

    # ------------------------------------------------------------------
    # 6. Ensamblar DataFrame de resultados
    # ------------------------------------------------------------------
    out = pd.DataFrame({
        "fecha": df["fecha"].values,
        "doy": df["doy"].values,
        "Tmin": df["Tmin"].values,
        "Tmax": df["Tmax"].values,
        "RHmin": df["RHmin"].values,
        "RHmax": df["RHmax"].values,
        "P_mm": df["P"].values,
        "Rs": df["Rs"].values,
        "u2_ms": df["u2"].values / 3.6,
        "ETo": ETo,
        "riego_mm": riego,
        "Kcb": Kcb,
        "Ke": Ke,
        "Kc": Kc,
        "Ks": Ks,
        "ETc_estandar_mm_d": ETc_std,
        "ETc_no_estandar_mm_d": ETc_ns,
        "De_mm": De,
        "Dr_mm": Dr,
    })
    out["ETc_estandar_cum_mm"] = out["ETc_estandar_mm_d"].cumsum()
    out["ETc_no_estandar_cum_mm"] = out["ETc_no_estandar_mm_d"].cumsum()
    out["riego_cum_mm"] = out["riego_mm"].cumsum()
    return out


# ===========================================================================
# GRAFICADO
# ===========================================================================
def plot_results(res, out_png):
    fig, axs = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle("ETc por FAO-56 - Coeficiente doble Kc + riego al 50 % TAW\n"
                 "(Escenario único: reposición a FC, 15-09-2025 a 07-04-2026)",
                 fontsize=13, fontweight="bold")

    x = res["fecha"]

    # 1) ETo y ETc
    ax = axs[0, 0]
    ax.plot(x, res["ETo"], label="ETo (Penman-Monteith)", color="gray")
    ax.plot(x, res["ETc_estandar_mm_d"], label="ETc estándar (Kcb+Ke)·ETo",
            color="tab:blue")
    ax.plot(x, res["ETc_no_estandar_mm_d"], label="ETc no estándar (Ks·Kcb+Ke)·ETo",
            color="tab:red", linestyle="--")
    ax.set_ylabel("[mm/día]"); ax.set_title("Evapotranspiración diaria")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.tick_params(axis='x', rotation=45)

    # 2) Coeficientes Kc, Kcb, Ke, Ks
    ax = axs[0, 1]
    ax.plot(x, res["Kc"], label="Kc = Kcb + Ke", color="tab:green")
    ax.plot(x, res["Kcb"], label="Kcb (basal)", color="tab:olive")
    ax.plot(x, res["Ke"], label="Ke (suelo)", color="tab:cyan")
    ax.plot(x, res["Ks"], label="Ks (estrés hídrico)", color="tab:red", linestyle=":")
    ax.set_ylabel("[-]"); ax.set_title("Coeficientes de cultivo y estrés")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.tick_params(axis='x', rotation=45)

    # 3) Cumuladas
    ax = axs[1, 0]
    ax.plot(x, res["ETc_estandar_cum_mm"], label="ETc estándar acumulada",
            color="tab:blue")
    ax.plot(x, res["ETc_no_estandar_cum_mm"], label="ETc no estándar acumulada",
            color="tab:red", linestyle="--")
    ax.set_ylabel("[mm]"); ax.set_title("ETc acumulada en el ciclo")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.tick_params(axis='x', rotation=45)

    # 4) Precipitación + riego vs ETc no estándar
    ax = axs[1, 1]
    ax.bar(x, res["P_mm"], width=1.0, alpha=0.4, label="Precipitación [mm]", color="tab:blue")
    riego_dias = res["riego_mm"] > 0
    if riego_dias.any():
        ax.bar(x[riego_dias], res.loc[riego_dias, "riego_mm"], width=1.0,
               alpha=0.8, label="Riego [mm]", color="tab:orange")
    ax.plot(x, res["ETc_no_estandar_mm_d"], label="ETc no estándar [mm/día]",
            color="tab:red")
    ax.set_ylabel("[mm]"); ax.set_title("Precipitación, riego y ETc no estándar")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.tick_params(axis='x', rotation=45)

    # 5) Depleción del suelo y umbral de riego
    ax = axs[2, 0]
    ax.plot(x, res["De_mm"], label="De - capa superficial [mm]", color="tab:orange")
    ax.plot(x, res["Dr_mm"], label="Dr - zona radicular [mm]", color="tab:purple")
    TAW = 1000.0 * (0.28 - 0.13) * 0.80   # mismo valor que en CropSoil por defecto
    ax.axhline(0.5 * TAW, color="black", linestyle="--", linewidth=1,
               label="Umbral riego (50 % TAW)")
    ax.set_ylabel("[mm]"); ax.set_title("Depleción de agua en el suelo")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.tick_params(axis='x', rotation=45)

    # 6) Riego acumulado y Ks
    ax = axs[2, 1]
    ax.plot(x, res["riego_cum_mm"], color="tab:orange", label="Riego acumulado [mm]")
    ax.set_ylabel("[mm]"); ax.set_title("Riego acumulado aplicado")
    ax2 = ax.twinx()
    ax2.plot(x, res["Ks"], color="tab:red", label="Ks", linestyle=":")
    ax2.set_ylabel("[-]"); ax2.set_ylim(0, 1.05)
    ax2.axhline(1.0, color="black", linestyle="--", linewidth=0.6)
    ax.legend(fontsize=8, loc="upper left"); ax2.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3); ax.tick_params(axis='x', rotation=45)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Gráfico guardado en: {out_png}")


# ===========================================================================
# GRAFICO ADICIONAL 1: Kc y Kcb a lo largo de la temporada
# ===========================================================================
def plot_kc_season(res, out_png):
    """Evolución temporal de Kc (doble) y sus componentes Kcb y Ke."""
    fig, ax = plt.subplots(figsize=(14, 5))
    x = res["fecha"]

    ax.plot(x, res["Kc"], label="Kc = Kcb + Ke (doble)", color="tab:green", linewidth=2)
    ax.plot(x, res["Kcb"], label="Kcb (coeficiente basal)", color="tab:olive", linewidth=2)
    ax.plot(x, res["Ke"], label="Ke (evaporación del suelo)", color="tab:cyan", linewidth=2)
    ax.plot(x, res["Ks"], label="Ks (estrés hídrico)", color="tab:red",
            linestyle=":", linewidth=2)

    ax.set_ylabel("Coeficiente [-]")
    ax.set_title("Coeficiente de cultivo Kc y componentes a lo largo de la temporada")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.25)
    ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Gráfico Kc temporada guardado en: {out_png}")


# ===========================================================================
# GRAFICO ADICIONAL 2: Contenido de agua del suelo, riegos y precipitaciones
# ===========================================================================
def plot_soil_water_riego(res, out_png, TAW=120.0):
    """
    Contenido de agua en la zona radicular (TAW - Dr) junto con precipitaciones
    y riegos aplicados en el tiempo.
    """
    fig, ax = plt.subplots(figsize=(14, 5))
    x = res["fecha"]

    # Contenido de agua del suelo (zona radicular) [mm]
    contenido_agua = TAW - res["Dr_mm"]
    ax.plot(x, contenido_agua, label="Agua en zona radicular (TAW - Dr)",
            color="tab:purple", linewidth=2)

    # Umbral de riego (50 % TAW) y capacidad de campo (TAW)
    ax.axhline(TAW, color="black", linestyle="-", linewidth=1, label="Capacidad de campo (FC)")
    ax.axhline(0.5 * TAW, color="orange", linestyle="--", linewidth=1.5, label="Umbral riego (50 % TAW)")

    # Barras de precipitación y riego
    ax.bar(x, res["P_mm"], width=1.0, alpha=0.35, label="Precipitación [mm]", color="tab:blue")
    riego_dias = res["riego_mm"] > 0
    if riego_dias.any():
        ax.bar(x[riego_dias], res.loc[riego_dias, "riego_mm"], width=1.0,
               alpha=0.85, label="Riego [mm]", color="tab:orange")

    ax.set_ylabel("Agua en suelo [mm]")
    ax.set_title("Contenido de agua del suelo, precipitaciones y riegos aplicados")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Gráfico agua del suelo y riego guardado en: {out_png}")


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    base = Path(__file__).resolve().parent
    xlsx = base / "Datos_clima_P2.xlsx"

    print(f"Leyendo datos climáticos de: {xlsx}")
    df = pd.read_excel(xlsx, sheet_name="ETc no estándar")

    # Estandarizar nombres de columnas y parsear fechas
    df = df.rename(columns={
        df.columns[0]: "fecha",
        df.columns[1]: "Tmin",
        df.columns[2]: "Tmax",
        df.columns[3]: "RHmin",
        df.columns[4]: "RHmax",
        df.columns[5]: "P",       # precipitación [mm]
        df.columns[6]: "Rs",      # radiación solar [MJ/m2/día]
        df.columns[7]: "u2",      # velocidad del viento [km/h] (a 2 m)
    })
    df["fecha"] = pd.to_datetime(df["fecha"], format="%d-%m-%Y")
    df["doy"] = df["fecha"].dt.dayofyear

    print(f"Registros leídos: {len(df)}  |  {df['fecha'].min().date()} a {df['fecha'].max().date()}")

    # Parámetros del cultivo/suelo
    cs = CropSoil()
    TAW = 1000.0 * (cs.th_FC - cs.th_WP) * cs.Zr
    TEW = 1000.0 * (cs.th_FC - cs.th_01) * cs.Ze

    # Ejecutar escenario ÚNICO con riego (reposición a FC al 50 % TAW)
    res = run_scenario(df, cs, lat=LAT, z=Z,
                       riego_umbral_frac=0.5, aplicar_riego=True)

    # Guardar resultados CSV
    out_csv = base / "resultados_ETc_escenario_riego.csv"
    res.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"Resultados guardados en: {out_csv}")

    # Guardar gráfico principal
    out_png = base / "fig_ETc_escenario_riego.png"
    plot_results(res, out_png)

    # Graficos adicionales solicitados
    out_kc = base / "fig_Kc_temporada.png"
    plot_kc_season(res, out_kc)

    TAW = 1000.0 * (cs.th_FC - cs.th_WP) * cs.Zr
    out_sw = base / "fig_agua_suelo_riego.png"
    plot_soil_water_riego(res, out_sw, TAW=TAW) 

    # ========================= RESUMEN =========================
    print("\n" + "=" * 70)
    print("RESUMEN DEL ESCENARIO ÚNICO (FAO-56, doble Kc, riego 50 % TAW)")
    print("=" * 70)
    print(f"Periodo simulado : {res['fecha'].min().date()}  a  {res['fecha'].max().date()}  ({len(res)} dias)")
    print(f"Cultivo/suelo    : Kcb_ini={cs.Kcb_ini}, Kcb_mid={cs.Kcb_mid}, Kcb_end={cs.Kcb_end}")
    print(f"                 : TAW={TAW:.0f} mm, TEW={TEW:.0f} mm, Umbral riego=50 % TAW = {0.5*TAW:.0f} mm")
    print(f"Riego            : reposición a capacidad de campo (FC)")
    print()
    print(f"ETo media        : {res['ETo'].mean():.2f} mm/día   (máx {res['ETo'].max():.2f})")
    print(f"ETc estándar total      : {res['ETc_estandar_cum_mm'].iloc[-1]:.1f} mm")
    print(f"ETc no estándar total   : {res['ETc_no_estandar_cum_mm'].iloc[-1]:.1f} mm")
    print(f"Riego total aplicado    : {res['riego_cum_mm'].iloc[-1]:.1f} mm")
    print(f"Número de riegos        : {(res['riego_mm'] > 0).sum()} eventos")
    print(f"Días con estrés (Ks < 0.95): {(res['Ks'] < 0.95).sum()} de {len(res)}")
    print(f"Ks mínimo               : {res['Ks'].min():.3f}")
    print(f"Kc medio                : {res['Kc'].mean():.2f}   (Ke medio {res['Ke'].mean():.2f})")
    print("=" * 70)


if __name__ == "__main__":
    main()
