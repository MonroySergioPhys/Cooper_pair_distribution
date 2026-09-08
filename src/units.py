"""
units.py
--------
Conversiones de unidades utilizadas en el proyecto
de distribución de pares de Cooper.
"""

from __future__ import annotations

import numpy as np


# ============================================================
# CONSTANTES
# ============================================================

# 1 Rydberg en eV
RY_TO_EV = 13.605693009

# 1 Rydberg en meV
RY_TO_MEV = RY_TO_EV * 1000.0

# 1 cm^-1 en eV
CM1_TO_EV = 1.0 / 8065.544005

# 1 cm^-1 en meV
CM1_TO_MEV = CM1_TO_EV * 1000.0

# 1 eV en meV
EV_TO_MEV = 1000.0

# 1 meV en eV
MEV_TO_EV = 1.0 / 1000.0


# ============================================================
# RYDBERG
# ============================================================

def rydberg_to_ev(value):
    """
    Convierte Rydberg a eV.
    """
    return np.asarray(value) * RY_TO_EV


def rydberg_to_mev(value):
    """
    Convierte Rydberg a meV.
    """
    return np.asarray(value) * RY_TO_MEV


# ============================================================
# cm^-1
# ============================================================

def cm1_to_ev(value):
    """
    Convierte cm^-1 a eV.
    """
    return np.asarray(value) * CM1_TO_EV


def cm1_to_mev(value):
    """
    Convierte cm^-1 a meV.
    """
    return np.asarray(value) * CM1_TO_MEV


# ============================================================
# eV
# ============================================================

def ev_to_mev(value):
    """
    Convierte eV a meV.
    """
    return np.asarray(value) * EV_TO_MEV


def mev_to_ev(value):
    """
    Convierte meV a eV.
    """
    return np.asarray(value) * MEV_TO_EV


# ============================================================
# ENERGÍA RELATIVA AL NIVEL DE FERMI
# ============================================================

def energy_relative_to_fermi(
    energy_ev,
    fermi_energy_ev
):
    """
    Calcula epsilon = E - EF.

    Parameters
    ----------
    energy_ev : array-like
        Energía absoluta en eV.

    fermi_energy_ev : float
        Energía de Fermi en eV.

    Returns
    -------
    np.ndarray
        Energía relativa al nivel de Fermi en eV.
    """

    return np.asarray(energy_ev) - fermi_energy_ev


def energy_relative_to_fermi_mev(
    energy_ev,
    fermi_energy_ev
):
    """
    Calcula epsilon = E - EF directamente en meV.
    """

    return energy_relative_to_fermi(
        energy_ev,
        fermi_energy_ev
    ) * EV_TO_MEV