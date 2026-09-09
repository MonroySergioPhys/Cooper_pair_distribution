"""
preprocessing.py
----------------
Procesamiento y preparación de los datos para el cálculo
de la distribución de pares de Cooper.

Convención interna:
    - energías fonónicas: meV
    - frecuencia fonónica auxiliar: THz
    - energía electrónica relativa a EF: meV

Principio fundamental:
    Los datos originales contenidos en los CSV representan la
    información física de entrada. No se generan artificialmente
    nuevos puntos del DOS para aumentar su resolución.

    Cuando una integral necesita evaluar el DOS entre dos puntos
    tabulados, se utiliza interpolación lineal sobre los datos
    originales.
"""

from __future__ import annotations

import numpy as np


# ============================================================
# UTILIDADES
# ============================================================

def remove_nan_and_inf(x, *ys):
    """
    Elimina posiciones donde x o cualquiera de los arrays y
    contengan NaN o infinito.
    """
    x = np.asarray(x, dtype=float)

    mask = np.isfinite(x)

    arrays = [x]

    for y in ys:
        y = np.asarray(y, dtype=float)
        mask &= np.isfinite(y)
        arrays.append(y)

    return tuple(array[mask] for array in arrays)


def _sort_and_unique(x, *ys):
    """
    Ordena x y elimina valores duplicados de x.

    Para un valor repetido se conserva la primera aparición.
    """
    order = np.argsort(x)

    x = x[order]
    ys = [y[order] for y in ys]

    x, unique_indices = np.unique(
        x,
        return_index=True,
    )

    ys = [
        y[unique_indices]
        for y in ys
    ]

    return (x, *ys)


# ============================================================
# ELIASHBERG
# ============================================================

def preprocess_eliashberg(
    frequency_mev,
    a2f_total,
    add_zero=True,
    add_cutoff=True,
):
    """
    Procesa la función de Eliashberg.

    Parámetros
    ----------
    frequency_mev : array_like
        Frecuencia/energía fonónica en meV.

    a2f_total : array_like
        Función alpha^2 F(omega).

    add_zero : bool
        Añade (0, 0) si no está presente.

    add_cutoff : bool
        Añade (omega_c, 0) si es necesario.

    Returns
    -------
    dict
        frequency_mev
        a2f_total
        omega_c_mev
    """

    frequency_mev = np.asarray(
        frequency_mev,
        dtype=float,
    )

    a2f_total = np.asarray(
        a2f_total,
        dtype=float,
    )

    # --------------------------------------------------------
    # Limpieza
    # --------------------------------------------------------

    frequency_mev, a2f_total = remove_nan_and_inf(
        frequency_mev,
        a2f_total,
    )

    # Solo frecuencias físicas no negativas
    mask = frequency_mev >= 0.0

    frequency_mev = frequency_mev[mask]
    a2f_total = a2f_total[mask]

    # Los pequeños valores negativos se consideran
    # artefactos numéricos.
    a2f_total = np.maximum(
        a2f_total,
        0.0,
    )

    frequency_mev, a2f_total = _sort_and_unique(
        frequency_mev,
        a2f_total,
    )

    if len(frequency_mev) == 0:
        raise ValueError(
            "No hay datos válidos de Eliashberg."
        )

    # --------------------------------------------------------
    # Determinación de omega_c
    # --------------------------------------------------------

    zero_indices = np.flatnonzero(
        a2f_total == 0.0
    )

    if len(zero_indices) == 0:
        raise ValueError(
            "No se encontró un punto con alpha^2F = 0. "
            "No es posible determinar omega_c automáticamente."
        )

    omega_c = float(
        frequency_mev[zero_indices[0]]
    )

    # Conservar la región hasta omega_c
    mask = frequency_mev <= omega_c

    frequency_mev = frequency_mev[mask]
    a2f_total = a2f_total[mask]

    if len(frequency_mev) == 0:
        raise ValueError(
            "Los datos de Eliashberg quedaron vacíos."
        )

    # --------------------------------------------------------
    # Condición alpha^2 F(0) = 0
    # --------------------------------------------------------

    if add_zero:

        if frequency_mev[0] > 0.0:

            frequency_mev = np.insert(
                frequency_mev,
                0,
                0.0,
            )

            a2f_total = np.insert(
                a2f_total,
                0,
                0.0,
            )

        else:

            frequency_mev[0] = 0.0
            a2f_total[0] = 0.0

    # --------------------------------------------------------
    # Condición alpha^2 F(omega_c) = 0
    # --------------------------------------------------------

    if add_cutoff:

        if frequency_mev[-1] < omega_c:

            frequency_mev = np.append(
                frequency_mev,
                omega_c,
            )

            a2f_total = np.append(
                a2f_total,
                0.0,
            )

        else:

            frequency_mev[-1] = omega_c
            a2f_total[-1] = 0.0

    return {
        "frequency_mev": frequency_mev,
        "a2f_total": a2f_total,
        "omega_c_mev": omega_c,
    }


# ============================================================
# PHONON DOS
# ============================================================

def preprocess_phonon_dos(
    frequency_mev,
    dos,
    omega_c_mev,
    pdos=None,
):
    """
    Procesa la densidad de estados fonónica.

    Los datos nuevos ya están expresados en meV.

    No se fuerza artificialmente:

        N_ph(omega_c) = 0

    porque el archivo puede no contener exactamente un punto
    en omega_c.
    """

    frequency_mev = np.asarray(
        frequency_mev,
        dtype=float,
    )

    dos = np.asarray(
        dos,
        dtype=float,
    )

    # --------------------------------------------------------
    # Limpieza
    # --------------------------------------------------------

    if pdos is not None:

        pdos = np.asarray(
            pdos,
            dtype=float,
        )

        frequency_mev, dos, pdos = remove_nan_and_inf(
            frequency_mev,
            dos,
            pdos,
        )

    else:

        frequency_mev, dos = remove_nan_and_inf(
            frequency_mev,
            dos,
        )

    # Solo frecuencias no negativas
    mask = frequency_mev >= 0.0

    frequency_mev = frequency_mev[mask]
    dos = dos[mask]

    if pdos is not None:
        pdos = pdos[mask]

    # El DOS no puede ser negativo físicamente
    dos = np.maximum(
        dos,
        0.0,
    )

    if pdos is not None:

        pdos = np.maximum(
            pdos,
            0.0,
        )

    # Ordenar y eliminar duplicados
    if pdos is not None:

        frequency_mev, dos, pdos = _sort_and_unique(
            frequency_mev,
            dos,
            pdos,
        )

    else:

        frequency_mev, dos = _sort_and_unique(
            frequency_mev,
            dos,
        )

    # --------------------------------------------------------
    # Región física
    # --------------------------------------------------------

    mask = frequency_mev <= float(
        omega_c_mev
    )

    frequency_mev = frequency_mev[mask]
    dos = dos[mask]

    if pdos is not None:
        pdos = pdos[mask]

    if len(frequency_mev) == 0:
        raise ValueError(
            "El DOS fonónico no contiene datos "
            "dentro del cutoff."
        )

    # --------------------------------------------------------
    # Condición N_ph(0) = 0
    # --------------------------------------------------------

    if frequency_mev[0] > 0.0:

        frequency_mev = np.insert(
            frequency_mev,
            0,
            0.0,
        )

        dos = np.insert(
            dos,
            0,
            0.0,
        )

        if pdos is not None:

            pdos = np.insert(
                pdos,
                0,
                0.0,
            )

    else:

        frequency_mev[0] = 0.0
        dos[0] = 0.0

        if pdos is not None:
            pdos[0] = 0.0

    result = {
        "frequency_mev": frequency_mev,
        "dos": dos,
        "omega_c_mev": float(omega_c_mev),
    }

    if pdos is not None:
        result["pdos"] = pdos

    return result


# ============================================================
# ELECTRONIC DOS
# ============================================================

def preprocess_electronic_dos(
    epsilon_mev,
    dos,
    omega_c_mev,
):
    """
    Procesa la densidad de estados electrónica.

    Se utiliza:

        epsilon = E - E_F

    y epsilon está expresado en meV.

    IMPORTANTE
    ----------
    NO se recortan los datos originales al intervalo
    [-omega_c, omega_c].

    Esto es necesario cuando la resolución del DOS electrónico
    es mayor que la ventana energética de interés.

    Por ejemplo, para Pb puede ocurrir:

        omega_c ~ 9.87 meV

    mientras que los datos originales contienen:

        -37 meV
        +13 meV

    Aunque ninguno de estos puntos está dentro de la ventana,
    ambos son necesarios para interpolar N_e(epsilon) dentro
    de dicha ventana.

    No se generan nuevos datos físicos. La interpolación se
    realiza exclusivamente entre puntos originalmente tabulados.
    """

    epsilon_mev = np.asarray(
        epsilon_mev,
        dtype=float,
    )

    dos = np.asarray(
        dos,
        dtype=float,
    )

    # --------------------------------------------------------
    # Limpieza
    # --------------------------------------------------------

    epsilon_mev, dos = remove_nan_and_inf(
        epsilon_mev,
        dos,
    )

    # El DOS no puede ser negativo físicamente.
    # Los valores negativos se interpretan como artefactos
    # numéricos del archivo.
    dos = np.maximum(
        dos,
        0.0,
    )

    # Ordenar y eliminar duplicados
    epsilon_mev, dos = _sort_and_unique(
        epsilon_mev,
        dos,
    )

    if len(epsilon_mev) == 0:
        raise ValueError(
            "No hay datos válidos de DOS electrónico."
        )

    omega_c_mev = float(
        omega_c_mev
    )

    # --------------------------------------------------------
    # Verificación de cobertura alrededor de EF
    # --------------------------------------------------------

    if epsilon_mev[0] > 0.0:

        raise ValueError(
            "El DOS electrónico no contiene datos por debajo "
            "de E_F. No es posible cubrir el intervalo físico."
        )

    if epsilon_mev[-1] < 0.0:

        raise ValueError(
            "El DOS electrónico no contiene datos por encima "
            "de E_F. No es posible cubrir el intervalo físico."
        )

    # --------------------------------------------------------
    # Interpolador
    # --------------------------------------------------------

    def dos_interpolator(epsilon):

        epsilon = np.asarray(
            epsilon,
            dtype=float,
        )

        return np.interp(
            epsilon,
            epsilon_mev,
            dos,
            left=0.0,
            right=0.0,
        )

    # --------------------------------------------------------
    # DOS en el nivel de Fermi
    # --------------------------------------------------------

    dos_at_fermi = float(
        dos_interpolator(0.0)
    )

    # --------------------------------------------------------
    # Resolución original del DOS
    # --------------------------------------------------------

    if len(epsilon_mev) > 1:

        delta_e_mev = np.diff(
            epsilon_mev
        )

        resolution_min_mev = float(
            np.min(delta_e_mev)
        )

        resolution_max_mev = float(
            np.max(delta_e_mev)
        )

    else:

        resolution_min_mev = np.nan
        resolution_max_mev = np.nan

    # --------------------------------------------------------
    # Puntos originales dentro de la ventana física
    # --------------------------------------------------------

    inside_window = (
        (epsilon_mev >= -omega_c_mev)
        &
        (epsilon_mev <= omega_c_mev)
    )

    n_points_inside = int(
        np.count_nonzero(
            inside_window
        )
    )

    # --------------------------------------------------------
    # Distancia a los puntos originales más cercanos
    # --------------------------------------------------------

    below = epsilon_mev[
        epsilon_mev <= 0.0
    ]

    above = epsilon_mev[
        epsilon_mev >= 0.0
    ]

    if len(below) > 0:
        epsilon_below_fermi = float(
            below[-1]
        )
    else:
        epsilon_below_fermi = np.nan

    if len(above) > 0:
        epsilon_above_fermi = float(
            above[0]
        )
    else:
        epsilon_above_fermi = np.nan

    # --------------------------------------------------------
    # Resultado
    # --------------------------------------------------------

    return {
        # ----------------------------------------------------
        # Datos originales
        # ----------------------------------------------------

        "epsilon_mev": epsilon_mev,
        "dos": dos,

        # ----------------------------------------------------
        # Parámetros físicos
        # ----------------------------------------------------

        "fermi_energy_ev": 0.0,

        "omega_c_mev": omega_c_mev,

        "epsilon_window_mev": (
            -omega_c_mev,
            omega_c_mev,
        ),

        # ----------------------------------------------------
        # Interpolación
        # ----------------------------------------------------

        "dos_interpolator": dos_interpolator,

        "dos_at_fermi": dos_at_fermi,

        # ----------------------------------------------------
        # Diagnóstico de resolución
        # ----------------------------------------------------

        "resolution_min_mev": (
            resolution_min_mev
        ),

        "resolution_max_mev": (
            resolution_max_mev
        ),

        "n_points_inside_window": (
            n_points_inside
        ),

        "has_direct_points_inside_window": (
            n_points_inside > 0
        ),

        "epsilon_below_fermi_mev": (
            epsilon_below_fermi
        ),

        "epsilon_above_fermi_mev": (
            epsilon_above_fermi
        ),
    }


# ============================================================
# PROCESAMIENTO COMPLETO
# ============================================================

def preprocess_aluminum(
    eliashberg_data,
    phonon_data,
    electronic_data,
):
    """
    Procesa los tres conjuntos de datos.

    El nombre de la función se conserva por compatibilidad
    con la estructura original del proyecto.
    """

    # ========================================================
    # ELIASHBERG
    # ========================================================

    eliashberg_frequency = eliashberg_data.get(
        "frequency_mev",
        eliashberg_data.get("frequency_ry"),
    )

    if eliashberg_frequency is None:

        raise KeyError(
            "Los datos de Eliashberg no contienen "
            "frequency_mev."
        )

    # Compatibilidad con datos antiguos
    if "frequency_mev" not in eliashberg_data:

        from .units import rydberg_to_mev

        eliashberg_frequency = rydberg_to_mev(
            eliashberg_frequency
        )

    eliashberg = preprocess_eliashberg(
        eliashberg_frequency,
        eliashberg_data["a2f_total"],
    )

    omega_c_mev = eliashberg[
        "omega_c_mev"
    ]

    # ========================================================
    # PHONON DOS
    # ========================================================

    phonon_frequency = phonon_data.get(
        "frequency_mev",
        phonon_data.get("frequency_cm1"),
    )

    if phonon_frequency is None:

        raise KeyError(
            "Los datos fonónicos no contienen "
            "frequency_mev."
        )

    # Compatibilidad con datos antiguos
    if "frequency_mev" not in phonon_data:

        from .units import cm1_to_mev

        phonon_frequency = cm1_to_mev(
            phonon_frequency
        )

    phonon = preprocess_phonon_dos(
        phonon_frequency,
        phonon_data["dos"],
        omega_c_mev,
        pdos=phonon_data.get("pdos"),
    )

    # ========================================================
    # ELECTRONIC DOS
    # ========================================================

    if "epsilon_mev" in electronic_data:

        epsilon_mev = electronic_data[
            "epsilon_mev"
        ]

    elif "energy_ev" in electronic_data:

        from .units import (
            energy_relative_to_fermi_mev
        )

        epsilon_mev = (
            energy_relative_to_fermi_mev(
                electronic_data["energy_ev"],
                electronic_data["fermi_energy_ev"],
            )
        )

    else:

        raise KeyError(
            "Los datos electrónicos no contienen "
            "epsilon_mev."
        )

    electronic = preprocess_electronic_dos(
        epsilon_mev,
        electronic_data["dos"],
        omega_c_mev,
    )

    # ========================================================
    # RESULTADO FINAL
    # ========================================================

    return {
        "eliashberg": eliashberg,

        "phonon": phonon,

        "electronic": electronic,

        "superconducting_parameters": {
            "omega_c_mev": omega_c_mev,

            "fermi_energy_ev": electronic_data.get(
                "fermi_energy_ev",
                0.0,
            ),

            "lambda": eliashberg_data.get(
                "lambda"
            ),

            "Delta": eliashberg_data.get(
                "Delta"
            ),
        },
    }