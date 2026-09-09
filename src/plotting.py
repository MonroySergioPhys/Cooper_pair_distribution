"""
plotting.py
-----------
Funciones para visualizar los datos procesados del proyecto
de distribución de pares de Cooper.

Convención de energía:
    meV
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def _prepare_xy(x, y):
    """Convierte a arrays y elimina NaN/inf."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)

    return x[mask], y[mask]


def plot_eliashberg(
    frequency_mev,
    a2f,
    omega_c_mev,
    material="Al",
    savepath=None,
):
    """Grafica la función de Eliashberg alpha^2 F(omega)."""

    frequency_mev, a2f = _prepare_xy(
        frequency_mev,
        a2f,
    )

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(
        frequency_mev,
        a2f,
        linewidth=2,
        label=r"$\alpha^2F(\omega)$",
    )

    ax.axvline(
        omega_c_mev,
        linestyle="--",
        linewidth=1.5,
        label=fr"$\omega_c = {omega_c_mev:.2f}$ meV",
    )

    ax.set_xlabel(
        r"$\omega$ (meV)",
        fontsize=12,
    )

    ax.set_ylabel(
        r"$\alpha^2F(\omega)$",
        fontsize=12,
    )

    ax.set_title(
        f"Función de Eliashberg — {material}",
        fontsize=14,
    )

    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(
            savepath,
            dpi=300,
            bbox_inches="tight",
        )

    return fig, ax


def plot_phonon_dos(
    frequency_mev,
    dos,
    omega_c_mev,
    material="Al",
    pdos=None,
    savepath=None,
):
    """
    Grafica la densidad de estados fonónica N_ph(omega).

    ``frequency_mev`` ya está en meV; no se realiza ninguna
    conversión durante la visualización.
    """

    frequency_mev, dos = _prepare_xy(
        frequency_mev,
        dos,
    )

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(
        frequency_mev,
        dos,
        linewidth=2,
        label=r"$N_{ph}(\omega)$",
    )

    # El PDOS queda disponible para futuras visualizaciones,
    # pero no se dibuja automáticamente para no mezclar magnitudes.
    if pdos is not None:
        pdos = np.asarray(pdos, dtype=float)

    ax.axvline(
        omega_c_mev,
        linestyle="--",
        linewidth=1.5,
        label=fr"$\omega_c = {omega_c_mev:.2f}$ meV",
    )

    ax.set_xlabel(
        r"$\omega$ (meV)",
        fontsize=12,
    )

    ax.set_ylabel(
        r"$N_{ph}(\omega)$",
        fontsize=12,
    )

    ax.set_title(
        f"Densidad de estados fonónica — {material}",
        fontsize=14,
    )

    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(
            savepath,
            dpi=300,
            bbox_inches="tight",
        )

    return fig, ax


def plot_electronic_dos(
    epsilon_mev,
    dos,
    omega_c_mev,
    material="Al",
    savepath=None,
):
    """
    Grafica la densidad de estados electrónica alrededor
    del nivel de Fermi.

    epsilon = E - E_F

    La energía ya se recibe directamente en meV.
    """

    epsilon_mev, dos = _prepare_xy(
        epsilon_mev,
        dos,
    )

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(
        epsilon_mev,
        dos,
        linewidth=2,
        label=r"$N_e(\epsilon)$",
    )

    ax.axvline(
        0.0,
        linestyle="--",
        linewidth=1.5,
        label=r"$E_F$",
    )

    ax.axvline(
        -omega_c_mev,
        linestyle=":",
        linewidth=1.5,
        label=r"$-\omega_c$",
    )

    ax.axvline(
        omega_c_mev,
        linestyle=":",
        linewidth=1.5,
        label=r"$+\omega_c$",
    )

    ax.set_xlim(
        -omega_c_mev,
        omega_c_mev,
    )

    ax.set_xlabel(
        r"$\epsilon = E-E_F$ (meV)",
        fontsize=12,
    )

    ax.set_ylabel(
        r"$N_e(\epsilon)$",
        fontsize=12,
    )

    ax.set_title(
        f"Densidad de estados electrónica — {material}",
        fontsize=14,
    )

    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(
            savepath,
            dpi=300,
            bbox_inches="tight",
        )

    return fig, ax
