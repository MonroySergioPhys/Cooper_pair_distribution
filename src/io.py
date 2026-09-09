"""
io.py
-----
Funciones para leer los archivos de datos utilizados en el
estudio de la distribución de pares de Cooper.

Nuevos CSV de Al:
    Al_alpha2F.csv
    Al_DOS_fonones.csv
    Al_DOS_electronica.csv

Convención para los nuevos datos:
    - Eliashberg: energia_meV
    - DOS fonónico: energia_meV
    - DOS electrónico: energia_meV = E - E_F

Los archivos legacy .dos siguen siendo compatibles.
"""

from __future__ import annotations

import csv
from pathlib import Path
import re
from typing import Optional

import numpy as np

from .units import (
    cm1_to_mev,
    rydberg_to_mev,
    thz_to_mev,
    ev_to_mev,
)


# ============================================================
# UTILIDADES GENERALES
# ============================================================

def _normalize_header_name(name: str) -> str:
    """Normaliza nombres de columnas para comparación robusta."""
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def _find_column(columns: list[str], aliases: list[str]) -> str | None:
    """Busca una columna por nombre normalizado."""
    normalized = {
        _normalize_header_name(column): column
        for column in columns
    }

    normalized_aliases = [
        _normalize_header_name(alias)
        for alias in aliases
    ]

    for alias in normalized_aliases:
        if alias in normalized:
            return normalized[alias]

    for column in columns:
        normalized_column = _normalize_header_name(column)

        if any(
            alias in normalized_column
            for alias in normalized_aliases
        ):
            return column

    return None


def _read_csv_rows(
    filepath: str | Path,
) -> tuple[list[str], list[dict[str, str]]]:
    """Lee un CSV con encabezados y devuelve sus filas."""
    filepath = Path(filepath)

    with filepath.open(
        "r",
        encoding="utf-8",
        errors="ignore",
        newline="",
    ) as handle:
        reader = csv.reader(handle)
        rows = [
            row
            for row in reader
            if row and any(cell.strip() for cell in row)
        ]

    if not rows:
        raise ValueError(
            f"El archivo CSV está vacío: {filepath}"
        )

    header = [cell.strip() for cell in rows[0]]
    records: list[dict[str, str]] = []

    for row in rows[1:]:
        record = {}

        for index, value in enumerate(row[:len(header)]):
            record[header[index]] = value.strip()

        if record:
            records.append(record)

    return header, records


def _read_lines(filepath: str | Path) -> list[str]:
    """Lee un archivo de texto y devuelve sus líneas."""
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(
            f"No se encontró el archivo: {filepath}"
        )

    if not filepath.is_file():
        raise ValueError(
            f"La ruta no corresponde a un archivo: {filepath}"
        )

    return filepath.read_text(
        encoding="utf-8",
        errors="ignore",
    ).splitlines()


def _numeric_tokens(line: str) -> list[float]:
    """Extrae números en notación decimal o científica."""
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

    matches = re.findall(
        pattern,
        line,
        flags=re.VERBOSE,
    )

    return [float(value) for value in matches]


def _estimate_lambda_from_a2f(
    frequency_mev: np.ndarray,
    a2f_total: np.ndarray,
) -> float:
    """
    Estima lambda mediante:

        lambda = 2 integral [alpha^2F(omega) / omega] domega

    usando omega en meV.
    """
    frequency_mev = np.asarray(
        frequency_mev,
        dtype=float,
    )
    a2f_total = np.asarray(
        a2f_total,
        dtype=float,
    )

    positive = frequency_mev > 0.0

    if not np.any(positive):
        return 0.0

    omega = frequency_mev[positive]
    alpha2f = a2f_total[positive]

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        ratio = alpha2f / omega

    return float(
        2.0 * np.trapezoid(ratio, omega)
    )


# ============================================================
# ELIASHBERG
# ============================================================

def read_eliashberg(filepath: str | Path) -> dict:
    """
    Lee la función de Eliashberg.

    Para los nuevos CSV se prioriza:

        energia_meV

    Si no existe, se acepta:

        frecuencia_THz

    y se convierte a meV.

    Para archivos legacy .dos, la frecuencia se interpreta
    originalmente en Ry y se convierte a meV.

    Returns
    -------
    dict
        {
            "frequency_mev": np.ndarray,
            "a2f_total": np.ndarray,
            "modes": np.ndarray,
            "lambda": float | None,
            "Delta": float | None
        }
    """

    filepath = Path(filepath)

    lambda_value: Optional[float] = None
    delta_value: Optional[float] = None

    # --------------------------------------------------------
    # CSV NUEVO
    # --------------------------------------------------------

    if filepath.suffix.lower() == ".csv":
        header, rows = _read_csv_rows(filepath)

        if not rows:
            raise ValueError(
                f"No se encontraron filas con datos en {filepath}"
            )

        energy_name = _find_column(
            header,
            [
                "energia_meV",
                "energy_meV",
                "energia_mev",
                "energy_mev",
            ],
        )

        thz_name = _find_column(
            header,
            [
                "frecuencia_THz",
                "frequency_THz",
                "frecuencia_thz",
                "frequency_thz",
            ],
        )

        a2f_name = _find_column(
            header,
            [
                "alpha2F",
                "alpha_2F",
                "a2f_total",
                "a2f",
                "alpha2f_total",
            ],
        )

        if a2f_name is None:
            raise ValueError(
                "No se encontró la columna alpha2F "
                f"en {filepath}"
            )

        # Preferimos energia_meV porque ya es la unidad interna.
        if energy_name is None and thz_name is None:
            raise ValueError(
                "El CSV de Eliashberg debe contener "
                "energia_meV o frecuencia_THz."
            )

        frequencies_mev = []
        a2f_total = []

        for row in rows:
            try:
                if energy_name is not None:
                    frequency_value = float(
                        row[energy_name]
                    )
                else:
                    frequency_value = float(
                        row[thz_name]
                    )
                    frequency_value = float(
                        thz_to_mev(frequency_value)
                    )

                a2f_value = float(row[a2f_name])

            except (
                KeyError,
                TypeError,
                ValueError,
            ):
                continue

            frequencies_mev.append(frequency_value)
            a2f_total.append(a2f_value)

        frequencies_mev = np.asarray(
            frequencies_mev,
            dtype=float,
        )
        a2f_total = np.asarray(
            a2f_total,
            dtype=float,
        )

        if len(frequencies_mev) == 0:
            raise ValueError(
                f"No se encontraron datos numéricos "
                f"de Eliashberg en {filepath}"
            )

        lambda_value = _estimate_lambda_from_a2f(
            frequencies_mev,
            a2f_total,
        )

        return {
            "frequency_mev": frequencies_mev,
            "a2f_total": a2f_total,
            "modes": np.empty((0, 0)),
            "lambda": lambda_value,
            "Delta": delta_value,
        }

    # --------------------------------------------------------
    # LEGACY .dos
    # --------------------------------------------------------

    lines = _read_lines(filepath)

    frequencies_ry = []
    a2f_total = []
    modes = []

    for line in lines:
        stripped = line.strip()

        if "lambda" in stripped.lower():
            match_lambda = re.search(
                r"lambda\s*=\s*([-+0-9.eE]+)",
                stripped,
                flags=re.IGNORECASE,
            )

            match_delta = re.search(
                r"Delta\s*=\s*([-+0-9.eE]+)",
                stripped,
                flags=re.IGNORECASE,
            )

            if match_lambda:
                lambda_value = float(
                    match_lambda.group(1)
                )

            if match_delta:
                delta_value = float(
                    match_delta.group(1)
                )

            continue

        if not stripped or stripped.startswith("#"):
            continue

        values = _numeric_tokens(stripped)

        if len(values) < 2:
            continue

        frequencies_ry.append(values[0])
        a2f_total.append(values[1])

        if len(values) > 2:
            modes.append(values[2:])

    frequencies_ry = np.asarray(
        frequencies_ry,
        dtype=float,
    )
    a2f_total = np.asarray(
        a2f_total,
        dtype=float,
    )

    if modes:
        max_modes = max(
            len(row)
            for row in modes
        )

        modes_array = np.full(
            (len(modes), max_modes),
            np.nan,
            dtype=float,
        )

        for i, row in enumerate(modes):
            modes_array[
                i,
                :len(row)
            ] = row
    else:
        modes_array = np.empty((0, 0))

    if len(frequencies_ry) == 0:
        raise ValueError(
            f"No se encontraron datos numéricos "
            f"de Eliashberg en {filepath}"
        )

    frequencies_mev = rydberg_to_mev(
        frequencies_ry
    )

    if lambda_value is None:
        lambda_value = _estimate_lambda_from_a2f(
            frequencies_mev,
            a2f_total,
        )

    return {
        "frequency_mev": frequencies_mev,
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
    Lee la densidad de estados fonónica.

    Para los nuevos CSV se prioriza:

        energia_meV

    Si no existe, se acepta:

        frecuencia_THz

    Para archivos legacy .dos se utiliza frecuencia en cm^-1.

    Returns
    -------
    dict
        {
            "frequency_mev": np.ndarray,
            "dos": np.ndarray,
            "pdos": np.ndarray | None
        }
    """

    filepath = Path(filepath)

    # --------------------------------------------------------
    # CSV NUEVO
    # --------------------------------------------------------

    if filepath.suffix.lower() == ".csv":
        header, rows = _read_csv_rows(filepath)

        if not rows:
            raise ValueError(
                f"No se encontraron filas con datos en {filepath}"
            )

        energy_name = _find_column(
            header,
            [
                "energia_meV",
                "energy_meV",
                "energia_mev",
                "energy_mev",
            ],
        )

        thz_name = _find_column(
            header,
            [
                "frecuencia_THz",
                "frequency_THz",
                "frecuencia_thz",
                "frequency_thz",
            ],
        )

        dos_name = _find_column(
            header,
            [
                "DOS_fonones",
                "dos_fonones",
                "dos_fononico",
                "phonon_dos",
                "ph_dos",
                "dos",
            ],
        )

        pdos_name = _find_column(
            header,
            [
                "PDOS",
                "pdos",
                "pdos_fonones",
            ],
        )

        if dos_name is None:
            raise ValueError(
                "No se encontró la columna DOS_fonones "
                f"en {filepath}"
            )

        if energy_name is None and thz_name is None:
            raise ValueError(
                "El CSV fonónico debe contener "
                "energia_meV o frecuencia_THz."
            )

        frequency_mev = []
        dos = []
        pdos = []

        for row in rows:
            try:
                if energy_name is not None:
                    frequency_value = float(
                        row[energy_name]
                    )
                else:
                    frequency_value = float(
                        row[thz_name]
                    )
                    frequency_value = float(
                        thz_to_mev(frequency_value)
                    )

                dos_value = float(row[dos_name])

                if (
                    pdos_name is not None
                    and pdos_name in row
                ):
                    pdos_value = float(
                        row[pdos_name]
                    )
                else:
                    pdos_value = np.nan

            except (
                KeyError,
                TypeError,
                ValueError,
            ):
                continue

            frequency_mev.append(frequency_value)
            dos.append(dos_value)
            pdos.append(pdos_value)

        frequency_mev = np.asarray(
            frequency_mev,
            dtype=float,
        )
        dos = np.asarray(
            dos,
            dtype=float,
        )
        pdos = np.asarray(
            pdos,
            dtype=float,
        )

        if len(frequency_mev) == 0:
            raise ValueError(
                f"No se encontraron datos numéricos "
                f"de DOS fonónico en {filepath}"
            )

        if not np.any(np.isfinite(pdos)):
            pdos = None

        return {
            "frequency_mev": frequency_mev,
            "dos": dos,
            "pdos": pdos,
        }

    # --------------------------------------------------------
    # LEGACY .dos
    # --------------------------------------------------------

    lines = _read_lines(filepath)

    frequency_cm1 = []
    dos = []
    pdos = []

    for line in lines:
        stripped = line.strip()

        if not stripped or stripped.startswith("#"):
            continue

        values = _numeric_tokens(stripped)

        if len(values) < 2:
            continue

        frequency_cm1.append(values[0])
        dos.append(values[1])

        if len(values) >= 3:
            pdos.append(values[2])
        else:
            pdos.append(np.nan)

    frequency_cm1 = np.asarray(
        frequency_cm1,
        dtype=float,
    )
    dos = np.asarray(
        dos,
        dtype=float,
    )
    pdos = np.asarray(
        pdos,
        dtype=float,
    )

    if len(frequency_cm1) == 0:
        raise ValueError(
            f"No se encontraron datos numéricos "
            f"de DOS fonónico en {filepath}"
        )

    if not np.any(np.isfinite(pdos)):
        pdos = None

    return {
        "frequency_mev": cm1_to_mev(
            frequency_cm1
        ),
        "dos": dos,
        "pdos": pdos,
    }


# ============================================================
# ELECTRONIC DOS
# ============================================================

def read_electronic_dos(filepath: str | Path) -> dict:
    """
    Lee el DOS electrónico.

    Para los nuevos CSV se utiliza directamente:

        energia_meV

    donde:

        energia_meV = E - E_F

    También se reconoce E_minus_EF_eV como alternativa y se
    convierte a meV.

    En consecuencia, para los nuevos datos:

        E_F = 0

    dentro del sistema de coordenadas relativo al nivel de Fermi.

    Returns
    -------
    dict
        {
            "epsilon_mev": np.ndarray,
            "dos": np.ndarray,
            "integrated_dos": np.ndarray,
            "fermi_energy_ev": float
        }
    """

    filepath = Path(filepath)

    # --------------------------------------------------------
    # CSV NUEVO
    # --------------------------------------------------------

    if filepath.suffix.lower() == ".csv":
        header, rows = _read_csv_rows(filepath)

        if not rows:
            raise ValueError(
                f"No se encontraron filas con datos en {filepath}"
            )

        energy_mev_name = _find_column(
            header,
            [
                "energia_meV",
                "energy_meV",
                "energia_mev",
                "energy_mev",
            ],
        )

        energy_ev_name = _find_column(
            header,
            [
                "E_minus_EF_eV",
                "E_minus_EF",
                "epsilon_eV",
                "energy_eV",
            ],
        )

        dos_name = _find_column(
            header,
            [
                "DOS_electronica",
                "dos_electronica",
                "electronic_dos",
                "n_e",
                "dos",
            ],
        )

        integrated_name = _find_column(
            header,
            [
                "integrated_dos",
                "int_dos",
                "integrated",
                "int",
            ],
        )

        if dos_name is None:
            raise ValueError(
                "No se encontró la columna DOS_electronica "
                f"en {filepath}"
            )

        if (
            energy_mev_name is None
            and energy_ev_name is None
        ):
            raise ValueError(
                "El CSV electrónico debe contener "
                "energia_meV o E_minus_EF_eV."
            )

        epsilon_mev = []
        dos = []
        integrated_dos = []

        for row in rows:
            try:
                if energy_mev_name is not None:
                    epsilon_value = float(
                        row[energy_mev_name]
                    )
                else:
                    epsilon_value = float(
                        ev_to_mev(
                            float(row[energy_ev_name])
                        )
                    )

                dos_value = float(
                    row[dos_name]
                )

                if (
                    integrated_name is not None
                    and integrated_name in row
                ):
                    integrated_value = float(
                        row[integrated_name]
                    )
                else:
                    integrated_value = np.nan

            except (
                KeyError,
                TypeError,
                ValueError,
            ):
                continue

            epsilon_mev.append(epsilon_value)
            dos.append(dos_value)
            integrated_dos.append(
                integrated_value
            )

        epsilon_mev = np.asarray(
            epsilon_mev,
            dtype=float,
        )
        dos = np.asarray(
            dos,
            dtype=float,
        )
        integrated_dos = np.asarray(
            integrated_dos,
            dtype=float,
        )

        if len(epsilon_mev) == 0:
            raise ValueError(
                f"No se encontraron datos numéricos "
                f"de DOS electrónico en {filepath}"
            )

        return {
            "epsilon_mev": epsilon_mev,
            "dos": dos,
            "integrated_dos": integrated_dos,
            "fermi_energy_ev": 0.0,
        }

    # --------------------------------------------------------
    # LEGACY .dos
    # --------------------------------------------------------

    lines = _read_lines(filepath)

    energy_ev = []
    dos = []
    integrated_dos = []

    fermi_energy = None

    for line in lines:
        stripped = line.strip()

        if (
            "EFermi" in stripped
            or "Efermi" in stripped
        ):
            match = re.search(
                r"EFermi\s*=\s*([-+0-9.eE]+)\s*eV",
                stripped,
                flags=re.IGNORECASE,
            )

            if match:
                fermi_energy = float(
                    match.group(1)
                )

        if not stripped or stripped.startswith("#"):
            continue

        values = _numeric_tokens(stripped)

        if len(values) < 2:
            continue

        energy_ev.append(values[0])
        dos.append(values[1])

        if len(values) >= 3:
            integrated_dos.append(values[2])
        else:
            integrated_dos.append(np.nan)

    energy_ev = np.asarray(
        energy_ev,
        dtype=float,
    )
    dos = np.asarray(
        dos,
        dtype=float,
    )
    integrated_dos = np.asarray(
        integrated_dos,
        dtype=float,
    )

    if len(energy_ev) == 0:
        raise ValueError(
            f"No se encontraron datos numéricos "
            f"de DOS electrónico en {filepath}"
        )

    if fermi_energy is None:
        raise ValueError(
            f"No se pudo encontrar EFermi en {filepath}"
        )

    epsilon_mev = ev_to_mev(
        energy_ev - fermi_energy
    )

    return {
        "epsilon_mev": epsilon_mev,
        "dos": dos,
        "integrated_dos": integrated_dos,
        "fermi_energy_ev": fermi_energy,
    }
