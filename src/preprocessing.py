"""
preprocessing.py
----------------
Procesamiento y preparación de los datos de Quantum ESPRESSO
para el cálculo de la distribución de pares de Cooper.

Material actual: Aluminio (Al)
"""

from __future__ import annotations

import numpy as np

from .units import (
    rydberg_to_mev,
    cm1_to_mev,
    energy_relative_to_fermi_mev,
)


# ============================================================
# UTILIDADES
# ============================================================

def remove_nan_and_inf(x, *ys):
    """
    Elimina posiciones donde x o cualquiera de los arrays y
    contengan NaN o infinito.

    Returns
    -------
    tuple
        (x_clean, y1_clean, y2_clean, ...)
    """

    x = np.asarray(x, dtype=float)

    mask = np.isfinite(x)

    for y in ys:
        y = np.asarray(y, dtype=float)
        mask &= np.isfinite(y)

    arrays = [x]

    for y in ys:
        arrays.append(np.asarray(y, dtype=float))

    return tuple(array[mask] for array in arrays)


# ============================================================
# ELIASHBERG
# ============================================================

def preprocess_eliashberg(
    frequency_ry,
    a2f_total,
    add_zero=True,
    add_cutoff=True,
):
    """
    Procesa la función de Eliashberg.

    Pasos:
        1. Conversión Ry -> meV.
        2. Eliminación de valores no finitos.
        3. Eliminación de frecuencias negativas.
        4. Ordenamiento por frecuencia.
        5. Eliminación de frecuencias duplicadas.
        6. Inclusión de (0, 0).
        7. Determinación de omega_c.
        8. Inclusión de (omega_c, 0).

    Parameters
    ----------
    frequency_ry : array-like
        Frecuencias en Rydberg.

    a2f_total : array-like
        alpha^2 F total.

    Returns
    -------
    dict
        Datos procesados.
    """

    frequency_ry = np.asarray(frequency_ry, dtype=float)
    a2f_total = np.asarray(a2f_total, dtype=float)

    # --------------------------------------------------------
    # Eliminar NaN / inf
    # --------------------------------------------------------

    frequency_ry, a2f_total = remove_nan_and_inf(
        frequency_ry,
        a2f_total
    )

    # --------------------------------------------------------
    # Convertir frecuencia a meV
    # --------------------------------------------------------

    frequency_mev = rydberg_to_mev(frequency_ry)

    # --------------------------------------------------------
    # Eliminar frecuencias negativas
    # --------------------------------------------------------

    mask = frequency_mev >= 0.0

    frequency_mev = frequency_mev[mask]
    a2f_total = a2f_total[mask]

    # --------------------------------------------------------
    # Eliminar valores negativos de alpha^2 F
    # --------------------------------------------------------

    a2f_total = np.maximum(a2f_total, 0.0)

    # --------------------------------------------------------
    # Ordenar
    # --------------------------------------------------------

    order = np.argsort(frequency_mev)

    frequency_mev = frequency_mev[order]
    a2f_total = a2f_total[order]

    # --------------------------------------------------------
    # Eliminar frecuencias duplicadas
    # --------------------------------------------------------

    frequency_mev, unique_indices = np.unique(
        frequency_mev,
        return_index=True
    )

    a2f_total = a2f_total[unique_indices]

    # --------------------------------------------------------
    # Determinar omega_c ORIGINAL
    #
    # El archivo de Al tiene un último punto con:
    #
    # alpha^2 F = 0
    #
    # y ese punto corresponde al cutoff.
    # --------------------------------------------------------

    zero_indices = np.where(a2f_total == 0.0)[0]

    if len(zero_indices) > 0:

        # Buscamos el primer cero después de haber
        # entrado en la región física.
        omega_c = frequency_mev[zero_indices[0]]

    else:

        raise ValueError(
            "No se encontró un punto con alpha^2F = 0. "
            "Se necesita determinar omega_c mediante "
            "extrapolación."
        )

    # --------------------------------------------------------
    # Mantener datos hasta omega_c
    # --------------------------------------------------------

    mask = frequency_mev <= omega_c

    frequency_mev = frequency_mev[mask]
    a2f_total = a2f_total[mask]

    # --------------------------------------------------------
    # Garantizar (0, 0)
    # --------------------------------------------------------

    if add_zero:

        if frequency_mev[0] > 0.0:

            frequency_mev = np.insert(
                frequency_mev,
                0,
                0.0
            )

            a2f_total = np.insert(
                a2f_total,
                0,
                0.0
            )

        else:

            frequency_mev[0] = 0.0
            a2f_total[0] = 0.0

    # --------------------------------------------------------
    # Garantizar (omega_c, 0)
    # --------------------------------------------------------

    if add_cutoff:

        if frequency_mev[-1] < omega_c:

            frequency_mev = np.append(
                frequency_mev,
                omega_c
            )

            a2f_total = np.append(
                a2f_total,
                0.0
            )

        else:

            frequency_mev[-1] = omega_c
            a2f_total[-1] = 0.0

    return {
        "frequency_mev": frequency_mev,
        "a2f_total": a2f_total,
        "omega_c_mev": float(omega_c),
    }


# ============================================================
# PHONON DOS
# ============================================================

def preprocess_phonon_dos(
    frequency_cm1,
    dos,
    omega_c_mev,
    pdos=None,
):
    """
    Procesa la densidad de estados fonónica.

    Pasos:
        1. cm^-1 -> meV.
        2. Eliminar valores no finitos.
        3. Eliminar frecuencias negativas.
        4. Eliminar DOS negativos.
        5. Ordenar.
        6. Recortar hasta omega_c.
        7. Añadir (0,0).
        8. Añadir (omega_c,0).

    El PDOS se conserva, pero no se utiliza como Nph principal.

    Returns
    -------
    dict
    """

    frequency_cm1 = np.asarray(
        frequency_cm1,
        dtype=float
    )

    dos = np.asarray(dos, dtype=float)

    if pdos is not None:
        pdos = np.asarray(pdos, dtype=float)

    # --------------------------------------------------------
    # Limpiar
    # --------------------------------------------------------

    if pdos is not None:

        frequency_cm1, dos, pdos = remove_nan_and_inf(
            frequency_cm1,
            dos,
            pdos
        )

    else:

        frequency_cm1, dos = remove_nan_and_inf(
            frequency_cm1,
            dos
        )

    # --------------------------------------------------------
    # Convertir a meV
    # --------------------------------------------------------

    frequency_mev = cm1_to_mev(frequency_cm1)

    # --------------------------------------------------------
    # Frecuencia >= 0
    # --------------------------------------------------------

    mask = frequency_mev >= 0.0

    frequency_mev = frequency_mev[mask]
    dos = dos[mask]

    if pdos is not None:
        pdos = pdos[mask]

    # --------------------------------------------------------
    # DOS no negativa
    # --------------------------------------------------------

    dos = np.maximum(dos, 0.0)

    # --------------------------------------------------------
    # Ordenar
    # --------------------------------------------------------

    order = np.argsort(frequency_mev)

    frequency_mev = frequency_mev[order]
    dos = dos[order]

    if pdos is not None:
        pdos = pdos[order]

    # --------------------------------------------------------
    # Eliminar duplicados
    # --------------------------------------------------------

    frequency_mev, unique_indices = np.unique(
        frequency_mev,
        return_index=True
    )

    dos = dos[unique_indices]

    if pdos is not None:
        pdos = pdos[unique_indices]

    # --------------------------------------------------------
    # Recortar hasta omega_c
    # --------------------------------------------------------

    mask = frequency_mev <= omega_c_mev

    frequency_mev = frequency_mev[mask]
    dos = dos[mask]

    if pdos is not None:
        pdos = pdos[mask]

    # --------------------------------------------------------
    # Añadir (0,0)
    # --------------------------------------------------------

    if len(frequency_mev) == 0 or frequency_mev[0] > 0.0:

        frequency_mev = np.insert(
            frequency_mev,
            0,
            0.0
        )

        dos = np.insert(
            dos,
            0,
            0.0
        )

        if pdos is not None:
            pdos = np.insert(
                pdos,
                0,
                0.0
            )

    else:

        frequency_mev[0] = 0.0
        dos[0] = 0.0

        if pdos is not None:
            pdos[0] = 0.0

    # --------------------------------------------------------
    # Añadir (omega_c, 0)
    # --------------------------------------------------------

    if frequency_mev[-1] < omega_c_mev:

        frequency_mev = np.append(
            frequency_mev,
            omega_c_mev
        )

        dos = np.append(
            dos,
            0.0
        )

        if pdos is not None:
            pdos = np.append(
                pdos,
                0.0
            )

    else:

        frequency_mev[-1] = omega_c_mev
        dos[-1] = 0.0

        if pdos is not None:
            pdos[-1] = 0.0

    result = {
        "frequency_mev": frequency_mev,
        "dos": dos,
    }

    if pdos is not None:
        result["pdos"] = pdos

    return result


# ============================================================
# ELECTRONIC DOS
# ============================================================

def preprocess_electronic_dos(
    energy_ev,
    dos,
    fermi_energy_ev,
    omega_c_mev,
):
    """
    Procesa el DOS electrónico.

    Se transforma:

        epsilon = E - EF

    y posteriormente se convierte a meV.

    Finalmente se conserva únicamente:

        -omega_c <= epsilon <= omega_c

    Parameters
    ----------
    energy_ev : array-like
        Energía absoluta E [eV].

    dos : array-like
        DOS electrónico.

    fermi_energy_ev : float
        Energía de Fermi [eV].

    omega_c_mev : float
        Cutoff fonónico [meV].

    Returns
    -------
    dict
    """

    energy_ev = np.asarray(
        energy_ev,
        dtype=float
    )

    dos = np.asarray(
        dos,
        dtype=float
    )

    # --------------------------------------------------------
    # Limpiar
    # --------------------------------------------------------

    energy_ev, dos = remove_nan_and_inf(
        energy_ev,
        dos
    )

    # --------------------------------------------------------
    # Epsilon = E - EF
    # --------------------------------------------------------

    epsilon_mev = energy_relative_to_fermi_mev(
        energy_ev,
        fermi_energy_ev
    )

    # --------------------------------------------------------
    # El DOS no debería ser negativo.
    #
    # El archivo puede contener pequeños valores negativos
    # numéricos. Los convertimos a cero.
    # --------------------------------------------------------

    dos = np.maximum(dos, 0.0)

    # --------------------------------------------------------
    # Recortar al intervalo:
    #
    # -omega_c <= epsilon <= omega_c
    # --------------------------------------------------------

    mask = (
        (epsilon_mev >= -omega_c_mev)
        &
        (epsilon_mev <= omega_c_mev)
    )

    epsilon_mev = epsilon_mev[mask]
    dos = dos[mask]

    # --------------------------------------------------------
    # Ordenar
    # --------------------------------------------------------

    order = np.argsort(epsilon_mev)

    epsilon_mev = epsilon_mev[order]
    dos = dos[order]

    # --------------------------------------------------------
    # Eliminar duplicados
    # --------------------------------------------------------

    epsilon_mev, unique_indices = np.unique(
        epsilon_mev,
        return_index=True
    )

    dos = dos[unique_indices]

    return {
        "epsilon_mev": epsilon_mev,
        "dos": dos,
        "fermi_energy_ev": float(fermi_energy_ev),
        "omega_c_mev": float(omega_c_mev),
    }


# ============================================================
# PROCESAMIENTO COMPLETO DE AL
# ============================================================

def preprocess_aluminum(
    eliashberg_data,
    phonon_data,
    electronic_data,
):
    """
    Procesa simultáneamente los tres conjuntos de datos
    correspondientes al aluminio.

    Parameters
    ----------
    eliashberg_data : dict
        Salida de read_eliashberg().

    phonon_data : dict
        Salida de read_phonon_dos().

    electronic_data : dict
        Salida de read_electronic_dos().

    Returns
    -------
    dict
        Diccionario con todos los datos procesados.
    """

    # --------------------------------------------------------
    # 1. Eliashberg
    # --------------------------------------------------------

    eliashberg = preprocess_eliashberg(
        eliashberg_data["frequency_ry"],
        eliashberg_data["a2f_total"],
    )

    omega_c_mev = eliashberg["omega_c_mev"]

    # --------------------------------------------------------
    # 2. Phonon DOS
    # --------------------------------------------------------

    phonon = preprocess_phonon_dos(
        phonon_data["frequency_cm1"],
        phonon_data["dos"],
        omega_c_mev,
        pdos=phonon_data.get("pdos"),
    )

    # --------------------------------------------------------
    # 3. Electronic DOS
    # --------------------------------------------------------

    electronic = preprocess_electronic_dos(
        electronic_data["energy_ev"],
        electronic_data["dos"],
        electronic_data["fermi_energy_ev"],
        omega_c_mev,
    )

    return {
        "eliashberg": eliashberg,
        "phonon": phonon,
        "electronic": electronic,

        "superconducting_parameters": {
            "omega_c_mev": omega_c_mev,
            "fermi_energy_ev": electronic_data["fermi_energy_ev"],
            "lambda": eliashberg_data.get("lambda"),
            "Delta": eliashberg_data.get("Delta"),
        },
    }