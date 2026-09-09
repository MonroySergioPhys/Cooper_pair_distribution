"""
units.py
--------
Conversiones de unidades utilizadas en el proyecto
de distribución de pares de Cooper.

Convención interna para los nuevos CSV:
    energía/frecuencia fonónica -> meV
    frecuencia auxiliar          -> THz
    energía electrónica relativa -> meV
"""

from __future__ import annotations

import numpy as np


# ============================================================
# CONSTANTES
# ============================================================

# Rydberg
RY_TO_EV = 13.605693009
RY_TO_MEV = RY_TO_EV * 1000.0

# cm^-1
CM1_TO_EV = 1.0 / 8065.544005
CM1_TO_MEV = CM1_TO_EV * 1000.0

# eV <-> meV
EV_TO_MEV = 1000.0
MEV_TO_EV = 1.0 / EV_TO_MEV

# THz <-> meV
# E = h f
# h = 4.135667696 meV THz^-1
THZ_TO_MEV = 4.135667696
MEV_TO_THZ = 1.0 / THZ_TO_MEV


# ============================================================
# RYDBERG
# ============================================================

def rydberg_to_ev(value):
    """Convierte Rydberg a eV."""
    return np.asarray(value, dtype=float) * RY_TO_EV


def rydberg_to_mev(value):
    """Convierte Rydberg a meV."""
    return np.asarray(value, dtype=float) * RY_TO_MEV


# ============================================================
# cm^-1
# ============================================================

def cm1_to_ev(value):
    """Convierte cm^-1 a eV."""
    return np.asarray(value, dtype=float) * CM1_TO_EV


def cm1_to_mev(value):
    """Convierte cm^-1 a meV."""
    return np.asarray(value, dtype=float) * CM1_TO_MEV


# ============================================================
# eV
# ============================================================

def ev_to_mev(value):
    """Convierte eV a meV."""
    return np.asarray(value, dtype=float) * EV_TO_MEV


def mev_to_ev(value):
    """Convierte meV a eV."""
    return np.asarray(value, dtype=float) * MEV_TO_EV


# ============================================================
# THz
# ============================================================

def thz_to_mev(value):
    """Convierte THz a meV mediante E = h f."""
    return np.asarray(value, dtype=float) * THZ_TO_MEV


def mev_to_thz(value):
    """Convierte meV a THz."""
    return np.asarray(value, dtype=float) * MEV_TO_THZ


# ============================================================
# ENERGÍA RELATIVA AL NIVEL DE FERMI
# ============================================================

def energy_relative_to_fermi(energy_ev, fermi_energy_ev):
    """
    Calcula epsilon = E - EF.

    Parámetros
    ----------
    energy_ev : array-like
        Energía absoluta en eV.
    fermi_energy_ev : float
        Energía de Fermi en eV.
    """
    return np.asarray(energy_ev, dtype=float) - float(fermi_energy_ev)


def energy_relative_to_fermi_mev(energy_ev, fermi_energy_ev):
    """
    Calcula epsilon = E - EF y devuelve el resultado en meV.
    """
    return ev_to_mev(
        energy_relative_to_fermi(energy_ev, fermi_energy_ev)
    )
