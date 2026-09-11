
"""
plotting.py
-----------
Funciones para generar las gráficas del proyecto
de distribución de pares de Cooper.

Este módulo contiene únicamente funciones relacionadas
con la visualización de los resultados.

La estética se encuentra centralizada para mantener
una apariencia consistente en todas las figuras.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. PALETA DE COLORES
# ============================================================

COLORS = {
    # Colores principales
    "blue": "#0077FF",        # Azul eléctrico
    "purple": "#7B2CBF",      # Morado
    "orange": "#FF6B00",      # Naranja
    "turquoise": "#00B8A9",   # Turquesa
    "yellow": "#FFB703",      # Amarillo
    "red": "#EF476F",         # Coral / rojo

    # Colores auxiliares
    "dark": "#202124",        # Texto y ejes
    "gray": "#6C757D",        # Líneas secundarias
}


# ============================================================
# 2. ESTILO GENERAL
# ============================================================

def _apply_style(ax):
    """
    Aplica el estilo general utilizado en todas las gráficas.
    """

    # --------------------------------------------------------
    # Fondo
    # --------------------------------------------------------

    ax.set_facecolor("white")

    # --------------------------------------------------------
    # Bordes
    # --------------------------------------------------------

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_color(COLORS["dark"])
    ax.spines["bottom"].set_color(COLORS["dark"])

    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    # --------------------------------------------------------
    # Ticks
    # --------------------------------------------------------

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=11,
        colors=COLORS["dark"],
        length=5,
        width=1.0,
    )

    # --------------------------------------------------------
    # Cuadrícula
    # --------------------------------------------------------

    ax.grid(
        True,
        which="major",
        linestyle="--",
        linewidth=0.7,
        alpha=0.20,
    )

    ax.set_axisbelow(True)


# ============================================================
# 3. PREPARACIÓN DE DATOS
# ============================================================

def _prepare_xy(x, y):
    """
    Convierte los datos a arrays numéricos y elimina
    valores que no sean finitos.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)

    return x[mask], y[mask]


# ============================================================
# 4. FINALIZACIÓN DE FIGURAS
# ============================================================

def _finish_figure(fig, save_path=None):
    """
    Ajusta la distribución de la figura y, opcionalmente,
    guarda la imagen.

    Devuelve dos valores para mantener compatibilidad con
    el código anterior del proyecto:

        fig, _ = plot_...
    """

    fig.tight_layout()

    if save_path is not None:

        save_path = Path(save_path)
        save_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fig.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )

    # Mantener compatibilidad con los notebooks
    return fig, None


# ============================================================
# 5. FUNCIÓN DE ELIASHBERG
# ============================================================

def plot_eliashberg(
    frequency_mev,
    a2f_total,
    omega_c_mev=None,
    material=None,
    save_path=None,
):
    """
    Grafica la función de Eliashberg:

        α²F(ω)

    Parámetros
    ----------
    frequency_mev : array-like
        Energía/frecuencia fonónica en meV.

    a2f_total : array-like
        Función de Eliashberg α²F(ω).

    omega_c_mev : float, optional
        Frecuencia de corte en meV.

    material : str, optional
        Nombre del material.

    save_path : str o Path, optional
        Ruta donde guardar la figura.
    """

    x, y = _prepare_xy(
        frequency_mev,
        a2f_total,
    )

    fig, ax = plt.subplots(
        figsize=(8.5, 5.5)
    )

    _apply_style(ax)

    # --------------------------------------------------------
    # Curva principal
    # --------------------------------------------------------

    ax.plot(
        x,
        y,
        color=COLORS["blue"],
        linewidth=2.8,
        solid_capstyle="round",
        label=r"$\alpha^2F(\omega)$",
    )

    # Área bajo la curva
    ax.fill_between(
        x,
        y,
        0,
        color=COLORS["blue"],
        alpha=0.10,
    )

    # --------------------------------------------------------
    # Frecuencia de corte
    # --------------------------------------------------------

    if omega_c_mev is not None:

        ax.axvline(
            omega_c_mev,
            color=COLORS["red"],
            linestyle="--",
            linewidth=2.0,
            label=rf"$\omega_c={omega_c_mev:.2f}$ meV",
        )

    # --------------------------------------------------------
    # Título
    # --------------------------------------------------------

    title = r"Función de Eliashberg $\alpha^2F(\omega)$"

    if material is not None:
        title += f" — {material}"

    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        color=COLORS["dark"],
        pad=12,
    )

    # --------------------------------------------------------
    # Ejes
    # --------------------------------------------------------

    ax.set_xlabel(
        r"Energía fonónica $\omega$ (meV)",
        fontsize=13,
        color=COLORS["dark"],
    )

    ax.set_ylabel(
        r"$\alpha^2F(\omega)$",
        fontsize=13,
        color=COLORS["dark"],
    )

    # --------------------------------------------------------
    # Leyenda
    # --------------------------------------------------------

    ax.legend(
        fontsize=10.5,
        frameon=False,
        loc="best",
    )

    return _finish_figure(
        fig,
        save_path,
    )


# ============================================================
# 6. DENSIDAD DE ESTADOS FONÓNICA
# ============================================================

def plot_phonon_dos(
    frequency_mev,
    dos,
    omega_c_mev=None,
    material=None,
    save_path=None,
):
    """
    Grafica la densidad de estados fonónica.
    """

    x, y = _prepare_xy(
        frequency_mev,
        dos,
    )

    fig, ax = plt.subplots(
        figsize=(8.5, 5.5)
    )

    _apply_style(ax)

    # --------------------------------------------------------
    # Curva
    # --------------------------------------------------------

    ax.plot(
        x,
        y,
        color=COLORS["purple"],
        linewidth=2.8,
        solid_capstyle="round",
        label="DOS fonónica",
    )

    # Área bajo la curva
    ax.fill_between(
        x,
        y,
        0,
        color=COLORS["purple"],
        alpha=0.10,
    )

    # --------------------------------------------------------
    # Frecuencia de corte
    # --------------------------------------------------------

    if omega_c_mev is not None:

        ax.axvline(
            omega_c_mev,
            color=COLORS["red"],
            linestyle="--",
            linewidth=2.0,
            label=rf"$\omega_c={omega_c_mev:.2f}$ meV",
        )

    # --------------------------------------------------------
    # Título
    # --------------------------------------------------------

    title = "Densidad de estados fonónica"

    if material is not None:
        title += f" — {material}"

    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        color=COLORS["dark"],
        pad=12,
    )

    # --------------------------------------------------------
    # Ejes
    # --------------------------------------------------------

    ax.set_xlabel(
        r"Energía fonónica $\omega$ (meV)",
        fontsize=13,
        color=COLORS["dark"],
    )

    ax.set_ylabel(
        "DOS fonónica",
        fontsize=13,
        color=COLORS["dark"],
    )

    # --------------------------------------------------------
    # Leyenda
    # --------------------------------------------------------

    ax.legend(
        fontsize=10.5,
        frameon=False,
        loc="best",
    )

    return _finish_figure(
        fig,
        save_path,
    )


# ============================================================
# 7. DENSIDAD DE ESTADOS ELECTRÓNICA
# ============================================================

def plot_electronic_dos(
    epsilon_mev,
    dos,
    omega_c_mev=None,
    material=None,
    save_path=None,
):
    """
    Grafica la densidad de estados electrónica alrededor
    del nivel de Fermi.
    """

    x, y = _prepare_xy(
        epsilon_mev,
        dos,
    )

    fig, ax = plt.subplots(
        figsize=(8.5, 5.5)
    )

    _apply_style(ax)

    # --------------------------------------------------------
    # Curva
    # --------------------------------------------------------

    ax.plot(
        x,
        y,
        color=COLORS["orange"],
        linewidth=2.8,
        solid_capstyle="round",
        label="DOS electrónica",
    )

    # Área bajo la curva
    ax.fill_between(
        x,
        y,
        0,
        color=COLORS["orange"],
        alpha=0.10,
    )

    # --------------------------------------------------------
    # Nivel de Fermi
    # --------------------------------------------------------

    ax.axvline(
        0,
        color=COLORS["turquoise"],
        linestyle="--",
        linewidth=2.0,
        label=r"$E_F$",
    )

    # --------------------------------------------------------
    # Ventana energética
    # --------------------------------------------------------

    if omega_c_mev is not None:

        ax.axvline(
            -omega_c_mev,
            color=COLORS["gray"],
            linestyle=":",
            linewidth=1.6,
        )

        ax.axvline(
            omega_c_mev,
            color=COLORS["gray"],
            linestyle=":",
            linewidth=1.6,
        )

        ax.set_xlim(
            -omega_c_mev,
            omega_c_mev,
        )

    # --------------------------------------------------------
    # Título
    # --------------------------------------------------------

    title = "Densidad de estados electrónica"

    if material is not None:
        title += f" — {material}"

    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        color=COLORS["dark"],
        pad=12,
    )

    # --------------------------------------------------------
    # Ejes
    # --------------------------------------------------------

    ax.set_xlabel(
        r"$\epsilon-E_F$ (meV)",
        fontsize=13,
        color=COLORS["dark"],
    )

    ax.set_ylabel(
        "DOS electrónica",
        fontsize=13,
        color=COLORS["dark"],
    )

    # --------------------------------------------------------
    # Leyenda
    # --------------------------------------------------------

    ax.legend(
        fontsize=10.5,
        frameon=False,
        loc="best",
    )

    return _finish_figure(
        fig,
        save_path,
    )


# ============================================================
# 8. DISTRIBUCIÓN DE PARES DE COOPER
# ============================================================

def plot_dcp(
    omega_mev,
    dcp,
    omega_c_mev=None,
    material=None,
    save_path=None,
):
    """
    Grafica la distribución de pares de Cooper:

        D_cp(ω)
    """

    x, y = _prepare_xy(
        omega_mev,
        dcp,
    )

    fig, ax = plt.subplots(
        figsize=(8.5, 5.5)
    )

    _apply_style(ax)

    # --------------------------------------------------------
    # Curva
    # --------------------------------------------------------

    ax.plot(
        x,
        y,
        color=COLORS["turquoise"],
        linewidth=2.8,
        solid_capstyle="round",
        label=r"$D_{cp}(\omega)$",
    )

    # Área
    ax.fill_between(
        x,
        y,
        0,
        color=COLORS["turquoise"],
        alpha=0.12,
    )

    # --------------------------------------------------------
    # Frecuencia de corte
    # --------------------------------------------------------

    if omega_c_mev is not None:

        ax.axvline(
            omega_c_mev,
            color=COLORS["red"],
            linestyle="--",
            linewidth=2.0,
            label=rf"$\omega_c={omega_c_mev:.2f}$ meV",
        )

    # --------------------------------------------------------
    # Título
    # --------------------------------------------------------

    title = (
        r"Distribución de pares de Cooper "
        r"$D_{cp}(\omega)$"
    )

    if material is not None:
        title += f" — {material}"

    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        color=COLORS["dark"],
        pad=12,
    )

    # --------------------------------------------------------
    # Ejes
    # --------------------------------------------------------

    ax.set_xlabel(
        r"Energía $\omega$ (meV)",
        fontsize=13,
        color=COLORS["dark"],
    )

    ax.set_ylabel(
        r"$D_{cp}(\omega)$",
        fontsize=13,
        color=COLORS["dark"],
    )

    # --------------------------------------------------------
    # Leyenda
    # --------------------------------------------------------

    ax.legend(
        fontsize=10.5,
        frameon=False,
        loc="best",
    )

    return _finish_figure(
        fig,
        save_path,
    )


# ============================================================
# 9. DISTRIBUCIÓN ENERGÉTICA eDcp
# ============================================================

def plot_edcp(
    epsilon_mev,
    edcp,
    material=None,
    save_path=None,
):
    """
    Grafica la distribución energética de pares de Cooper:

        eD_cp(ε)
    """

    x, y = _prepare_xy(
        epsilon_mev,
        edcp,
    )

    fig, ax = plt.subplots(
        figsize=(8.5, 5.5)
    )

    _apply_style(ax)

    # --------------------------------------------------------
    # Curva
    # --------------------------------------------------------

    ax.plot(
        x,
        y,
        color=COLORS["yellow"],
        linewidth=2.8,
        solid_capstyle="round",
        label=r"$eD_{cp}(\epsilon)$",
    )

    # Área
    ax.fill_between(
        x,
        y,
        0,
        color=COLORS["yellow"],
        alpha=0.12,
    )

    # --------------------------------------------------------
    # Nivel de Fermi
    # --------------------------------------------------------

    ax.axvline(
        0,
        color=COLORS["red"],
        linestyle="--",
        linewidth=2.0,
        label=r"$E_F$",
    )

    # --------------------------------------------------------
    # Título
    # --------------------------------------------------------

    title = (
        r"Distribución energética "
        r"$eD_{cp}(\epsilon)$"
    )

    if material is not None:
        title += f" — {material}"

    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        color=COLORS["dark"],
        pad=12,
    )

    # --------------------------------------------------------
    # Ejes
    # --------------------------------------------------------

    ax.set_xlabel(
        r"$\epsilon-E_F$ (meV)",
        fontsize=13,
        color=COLORS["dark"],
    )

    ax.set_ylabel(
        r"$eD_{cp}(\epsilon)$",
        fontsize=13,
        color=COLORS["dark"],
    )

    # --------------------------------------------------------
    # Leyenda
    # --------------------------------------------------------

    ax.legend(
        fontsize=10.5,
        frameon=False,
        loc="best",
    )

    return _finish_figure(
        fig,
        save_path,
    )
