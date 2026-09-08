"""
io.py
-----
Funciones para leer los archivos de salida de Quantum ESPRESSO
utilizados en el estudio de la distribución de pares de Cooper.

Material actual: Aluminio (Al)

Archivos:
    Ala2F.dos1 -> función de Eliashberg alpha^2 F(omega)
    Alph.dos   -> DOS fonónico
    Aldos.dos  -> DOS electrónico
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Optional

import numpy as np


# ============================================================
# UTILIDADES GENERALES
# ============================================================

def _read_lines(filepath: str | Path) -> list[str]:
    """
    Lee un archivo de texto y devuelve sus líneas.

    Parameters
    ----------
    filepath : str or Path
        Ruta al archivo.

    Returns
    -------
    list[str]
        Lista de líneas del archivo.
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"No se encontró el archivo: {filepath}")

    if not filepath.is_file():
        raise ValueError(f"La ruta no corresponde a un archivo: {filepath}")

    return filepath.read_text(encoding="utf-8", errors="ignore").splitlines()


def _numeric_tokens(line: str) -> list[float]:
    """
    Extrae todos los números escritos en notación decimal o científica
    de una línea.

    Ejemplos aceptados:
        0.364967E-05
        -1.6635790180E-06
        7.931
        -4.041
        1.0e+03
    """

    pattern = r"""
        [-+]?
        (?:
            (?:\d+\.\d*) |
            (?:\.\d+) |
            (?:\d+)
        )
        (?:
            [Ee][-+]?\d+
        )?
    """

    matches = re.findall(pattern, line, flags=re.VERBOSE)

    return [float(value) for value in matches]


# ============================================================
# ELIASHBERG
# ============================================================

def read_eliashberg(filepath: str | Path) -> dict:
    """
    Lee un archivo tipo XXa2F.dosN.

    Para Al:
        Ala2F.dos1

    Formato esperado:

        frequency[Ry]  a2F_total  a2F(mode1) ...

    Se utiliza la primera columna como frecuencia y la segunda
    como alpha^2 F total.

    También intenta recuperar:
        lambda
        Delta

    Returns
    -------
    dict
        {
            "frequency_ry": np.ndarray,
            "a2f_total": np.ndarray,
            "modes": np.ndarray,
            "lambda": float | None,
            "Delta": float | None
        }
    """

    lines = _read_lines(filepath)

    frequencies = []
    a2f_total = []
    modes = []

    lambda_value: Optional[float] = None
    delta_value: Optional[float] = None

    for line in lines:

        stripped = line.strip()

        # ----------------------------------------------------
        # Buscar lambda y Delta
        # ----------------------------------------------------

        if "lambda" in stripped.lower():
            match_lambda = re.search(
                r"lambda\s*=\s*([-+0-9.eE]+)",
                stripped,
                flags=re.IGNORECASE
            )

            match_delta = re.search(
                r"Delta\s*=\s*([-+0-9.eE]+)",
                stripped,
                flags=re.IGNORECASE
            )

            if match_lambda:
                lambda_value = float(match_lambda.group(1))

            if match_delta:
                delta_value = float(match_delta.group(1))

            continue

        # ----------------------------------------------------
        # Ignorar comentarios
        # ----------------------------------------------------

        if not stripped or stripped.startswith("#"):
            continue

        values = _numeric_tokens(stripped)

        # Necesitamos al menos frecuencia + a2F total
        if len(values) < 2:
            continue

        frequencies.append(values[0])
        a2f_total.append(values[1])

        # Guardamos los modos disponibles
        if len(values) > 2:
            modes.append(values[2:])

    frequencies = np.asarray(frequencies, dtype=float)
    a2f_total = np.asarray(a2f_total, dtype=float)

    # Convertimos los modos a array rectangular si existen
    if modes:
        max_modes = max(len(row) for row in modes)

        modes_array = np.full(
            (len(modes), max_modes),
            np.nan,
            dtype=float
        )

        for i, row in enumerate(modes):
            modes_array[i, :len(row)] = row
    else:
        modes_array = np.empty((0, 0))

    if len(frequencies) == 0:
        raise ValueError(
            f"No se encontraron datos numéricos de Eliashberg en {filepath}"
        )

    return {
        "frequency_ry": frequencies,
        "a2f_total": a2f_total,
        "modes": modes_array,
        "lambda": lambda_value,
        "Delta": delta_value,
    }


# ============================================================
# PHONON DOS
# ============================================================

def read_phonon_dos(filepath: str | Path) -> dict:
    """
    Lee un archivo tipo XXph.dos.

    Formato esperado:

        Frequency[cm^-1]   DOS   PDOS

    La estructura se interpreta como:

        columna 1 -> frecuencia [cm^-1]
        columna 2 -> DOS fonónico total
        columna 3 -> PDOS

    Returns
    -------
    dict
        {
            "frequency_cm1": np.ndarray,
            "dos": np.ndarray,
            "pdos": np.ndarray
        }
    """

    lines = _read_lines(filepath)

    frequency = []
    dos = []
    pdos = []

    for line in lines:

        stripped = line.strip()

        if not stripped or stripped.startswith("#"):
            continue

        values = _numeric_tokens(stripped)

        # Necesitamos frecuencia + DOS
        if len(values) < 2:
            continue

        frequency.append(values[0])
        dos.append(values[1])

        # PDOS es opcional
        if len(values) >= 3:
            pdos.append(values[2])
        else:
            pdos.append(np.nan)

    frequency = np.asarray(frequency, dtype=float)
    dos = np.asarray(dos, dtype=float)
    pdos = np.asarray(pdos, dtype=float)

    if len(frequency) == 0:
        raise ValueError(
            f"No se encontraron datos numéricos de DOS fonónico en {filepath}"
        )

    return {
        "frequency_cm1": frequency,
        "dos": dos,
        "pdos": pdos,
    }


# ============================================================
# ELECTRONIC DOS
# ============================================================

def read_electronic_dos(filepath: str | Path) -> dict:
    """
    Lee un archivo tipo XXdos.dos.

    Formato esperado:

        E(eV)   dos(E)   Int dos(E) EFermi = X eV

    Para Al:

        EFermi = 7.931 eV

    Se utilizan:

        columna 1 -> energía E [eV]
        columna 2 -> DOS electrónico

    La columna de DOS integrado se conserva.

    Returns
    -------
    dict
        {
            "energy_ev": np.ndarray,
            "dos": np.ndarray,
            "integrated_dos": np.ndarray,
            "fermi_energy_ev": float
        }
    """

    lines = _read_lines(filepath)

    energy = []
    dos = []
    integrated_dos = []

    fermi_energy = None

    for line in lines:

        stripped = line.strip()

        # ----------------------------------------------------
        # Buscar EFermi
        # ----------------------------------------------------

        if "EFermi" in stripped or "Efermi" in stripped:

            match = re.search(
                r"EFermi\s*=\s*([-+0-9.eE]+)\s*eV",
                stripped,
                flags=re.IGNORECASE
            )

            if match:
                fermi_energy = float(match.group(1))

        # ----------------------------------------------------
        # Ignorar comentarios
        # ----------------------------------------------------

        if not stripped or stripped.startswith("#"):
            continue

        values = _numeric_tokens(stripped)

        if len(values) < 2:
            continue

        energy.append(values[0])
        dos.append(values[1])

        if len(values) >= 3:
            integrated_dos.append(values[2])
        else:
            integrated_dos.append(np.nan)

    energy = np.asarray(energy, dtype=float)
    dos = np.asarray(dos, dtype=float)
    integrated_dos = np.asarray(integrated_dos, dtype=float)

    if len(energy) == 0:
        raise ValueError(
            f"No se encontraron datos numéricos de DOS electrónico en {filepath}"
        )

    if fermi_energy is None:
        raise ValueError(
            f"No se pudo encontrar EFermi en {filepath}"
        )

    return {
        "energy_ev": energy,
        "dos": dos,
        "integrated_dos": integrated_dos,
        "fermi_energy_ev": fermi_energy,
    }