# Cooper Pair Distribution Function

<p align="center">
  <img src="https://img.shields.io/badge/Physics-Computational%20Physics-145DA0?style=for-the-badge" alt="Computational Physics">
  <img src="https://img.shields.io/badge/Python-3.x-1E88E5?style=for-the-badge&logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/Superconductivity-Cooper%20Pairs-00A6A6?style=for-the-badge" alt="Superconductivity">
</p>

<p align="center">
  <b>Cálculo numérico de la función de distribución de pares de Cooper a partir de datos de Eliashberg y densidades de estados.</b>
</p>

---

## 📖 Descripción

Este proyecto implementa un procedimiento computacional para estudiar la **distribución de pares de Cooper** en materiales superconductores a partir de información electrónica y fonónica.

El cálculo se basa en la función de distribución de pares de Cooper \(D_{cp}(\omega,T_c)\) y en sus distribuciones energéticas asociadas \(eD_{cp}(\epsilon)\) y \(eD_{cp}(\epsilon')\).

El objetivo principal es analizar cómo las propiedades electrónicas y fonónicas de un material se relacionan con la formación de pares de Cooper y, particularmente, estudiar la relación entre la separación energética obtenida a partir de la distribución calculada y la brecha superconductora experimental \(\Delta(0)\).

Los materiales considerados actualmente son:

* **Al** — Aluminio
* **Hg** — Mercurio
* **Nb** — Niobio
* **Ta** — Tantalio
* **Pb** — Plomo

---

## 🎯 Objetivos

### Objetivo general

Implementar computacionalmente la función de distribución de pares de Cooper utilizando datos de la función de Eliashberg, la densidad de estados fonónica y la densidad de estados electrónica de diferentes materiales superconductores.

### Objetivos específicos

* Procesar y homogenizar los datos de entrada.
* Determinar numéricamente la frecuencia de corte \(\omega_c\).
* Incorporar la función de Eliashberg \(\alpha^2F(\omega)\).
* Procesar la densidad de estados fonónica \(N_{ph}(\omega)\).
* Procesar la densidad de estados electrónica \(N_e(\epsilon)\).
* Calcular la distribución de pares de Cooper \(D_{cp}(\omega,T_c)\).
* Calcular \(eD_{cp}(\epsilon)\) y \(eD_{cp}(\epsilon')\).
* Determinar las energías características \(\omega_\epsilon\) y \(\omega_{\epsilon'}\).
* Calcular

$$
\Delta_{Dcp}
=
\omega_\epsilon-\omega_{\epsilon'}.
$$

* Comparar \(\frac{1}{2}\Delta_{Dcp}\) con la brecha superconductora \(\Delta(0)\).

---

# 🧠 Fundamento teórico

## Función de Eliashberg

La interacción electrón-fonón se describe mediante la función espectral de Eliashberg

$$
\alpha^2F(\omega).
$$

A partir de ella se obtiene el acoplamiento electrón-fonón mediante

$$
\alpha^2(\omega)
=
\frac{\alpha^2F(\omega)}
{N_{ph}(\omega)}.
$$

---

## Distribuciones electrónicas

Para los electrones ocupados se utiliza

$$
g_e^o(\epsilon,T)
=
N_e(\epsilon)
\frac{1}
{e^{\beta(\epsilon-E_F)}+1},
$$

mientras que para los estados vacantes:

$$
g_e^v(\epsilon,T)
=
N_e(\epsilon)
\left(
1-
\frac{1}
{e^{\beta(\epsilon-E_F)}+1}
\right).
$$

---

## Distribuciones fonónicas

Para los fonones:

$$
g_p^n(\omega,T)
=
N_{ph}(\omega)
\frac{1}
{e^{\beta\omega}-1},
$$

y

$$
g_p^{n+1}(\omega,T)
=
N_{ph}(\omega)
\left(
1+
\frac{1}
{e^{\beta\omega}-1}
\right).
$$

En la implementación numérica, las energías fonónicas se expresan en **meV**.

---

## Distribución de pares de Cooper

La función principal calculada por el proyecto es

$$
D_{cp}(\omega,T_c)=
\int_{E_F-\omega_c}^{E_F+\omega_c}
\int_{E_F-\omega_c}^{E_F+\omega_c}
g_e^o(\epsilon,T_c)
g_e^v(\epsilon-\omega,T_c)
g_p^{n+1}(\omega,T_c)
$$

$$
\times
g_e^o(\epsilon',T_c)
g_e^v(\epsilon'+\omega,T_c)
g_p^n(\omega,T_c)
\alpha^2(\omega)
\,d\epsilon\,d\epsilon'.
$$

La cantidad integrada sobre las frecuencias fonónicas es

$$
N_{cp}
=
\int_0^{\omega_c}
D_{cp}(\omega,T_c)\,d\omega.
$$

---

## Distribuciones energéticas

También se calculan las distribuciones

$$
eD_{cp}(\epsilon)
$$

y

$$
eD_{cp}(\epsilon').
$$

Los máximos de estas funciones permiten obtener

$$
\omega_\epsilon
\qquad\text{y}\qquad
\omega_{\epsilon'}.
$$

A partir de estas energías se define

$$
\boxed{
\Delta_{Dcp}
=
\omega_\epsilon-\omega_{\epsilon'}
}
$$

y se compara

$$
\boxed{
\frac{1}{2}\Delta_{Dcp}
}
$$

con la brecha superconductora \(\Delta(0)\).

---

# 🧮 Metodología computacional

El procesamiento sigue el siguiente flujo:

```text
              Datos de entrada
                    │
        ┌───────────┼───────────┐
        │           │           │
        ▼           ▼           ▼
     α²F(ω)      DOS fonónica   DOS electrónica
        │           │           │
        ▼           ▼           ▼
   Preprocesamiento y conversión de unidades
                    │
                    ▼
              Determinación
               de ω_c
                    │
                    ▼
        ┌───────────────────────┐
        │ Distribuciones        │
        │ electrónicas/fonónicas│
        └───────────┬───────────┘
                    │
                    ▼
             D_cp(ω,T_c)
                    │
             ┌──────┴──────┐
             ▼             ▼
       eD_cp(ε)       eD_cp(ε')
             │             │
             ▼             ▼
          ω_ε          ω_ε'
             │             │
             └──────┬──────┘
                    ▼
              Δ_Dcp
                    │
                    ▼
           ½ Δ_Dcp vs Δ(0)
```

---

# 📂 Estructura del proyecto

```text
Cooper_pair_distribution/
│
├── data/
│   └── raw/
│       ├── Al/
│       │   ├── Al_alpha2F.csv
│       │   ├── Al_DOS_fonones.csv
│       │   └── Al_DOS_electronica.csv
│       │
│       ├── Hg/
│       ├── Nb/
│       ├── Ta/
│       └── Pb/
│
├── src/
│   ├── io.py
│   ├── units.py
│   ├── preprocessing.py
│   ├── plotting.py
│   └── cooper_pairs.py
│
├── notebooks/
│   ├── 01_exploracion_datos.ipynb
│   ├── 02_calculo_Dcp.ipynb
│   ├── 03_calculo_eDcp.ipynb
│   └── 04_resultados_analisis.ipynb
│
├── results/
│   ├── figures/
│   └── tables/
│
└── README.md
```

---

# 📊 Datos de entrada

Para cada material se utilizan tres conjuntos de datos.

### 1. Función de Eliashberg

```text
<MATERIAL>_alpha2F.csv
```

Contiene, entre otras variables:

```text
frecuencia_THz
energia_meV
alpha2F
```

La energía fonónica utilizada internamente está expresada en meV.

---

### 2. Densidad de estados fonónica

```text
<MATERIAL>_DOS_fonones.csv
```

Contiene:

```text
frecuencia_THz
energia_meV
DOS_fonones
```

Los valores negativos asociados a errores numéricos del procesamiento de los datos se eliminan mediante el preprocesamiento.

---

### 3. Densidad de estados electrónica

```text
<MATERIAL>_DOS_electronica.csv
```

Contiene la energía electrónica relativa al nivel de Fermi:

```text
E_minus_EF_eV
energia_meV
DOS_electronica
```

El nivel de Fermi se toma como

$$
E_F=0.
$$

---

# ⚙️ Unidades

La convención interna del proyecto es:

| Magnitud                     | Unidad |
| ---------------------------- | ------ |
| Energía electrónica          | meV    |
| Energía fonónica             | meV    |
| Frecuencia fonónica auxiliar | THz    |
| Temperatura                  | K      |

La conversión utilizada entre THz y meV es

$$
1\ \mathrm{THz}
\approx
4.135667696\ \mathrm{meV}.
$$

Las conversiones se encuentran centralizadas en:

```text
src/units.py
```

---

# 🔬 Preprocesamiento

Antes del cálculo se realizan diferentes operaciones:

* eliminación de valores no válidos;
* eliminación de frecuencias negativas;
* ordenamiento de los datos;
* eliminación de puntos duplicados;
* tratamiento de valores negativos espurios en las DOS;
* incorporación del origen cuando es necesario;
* determinación de \(\omega_c\);
* interpolación de la densidad de estados electrónica.

Un punto importante del procesamiento es que **no se generan datos físicos artificiales para aumentar la resolución de las DOS**.

Las mallas numéricas utilizadas para integrar las expresiones se distinguen de los puntos físicos proporcionados por los datos originales.

---

# 🚀 Ejecución

## 1. Clonar el repositorio

```bash
git clone https://github.com/MonroySergioPhys/Cooper_pair_distribution.git
cd Cooper_pair_distribution
```

---

## 2. Crear un entorno virtual

Se recomienda utilizar un entorno virtual de Python:

```bash
python -m venv .venv
```

### Linux/macOS

```bash
source .venv/bin/activate
```

### Windows

```powershell
.venv\Scripts\activate
```

---

## 3. Instalar dependencias

Instalar las principales dependencias científicas:

```bash
pip install numpy scipy pandas matplotlib jupyter
```

---

# 📓 Notebooks

El análisis está dividido en cuatro notebooks para mantener separado el procesamiento, el cálculo y la interpretación de resultados.

## `01_exploracion_datos.ipynb`

Permite:

* cargar los datos;
* verificar las unidades;
* inspeccionar las distribuciones;
* realizar el preprocesamiento;
* determinar \(\omega_c\);
* generar las primeras gráficas.

---

## `02_calculo_Dcp.ipynb`

Calcula

$$
D_{cp}(\omega,T_c)
$$

para cada material.

También obtiene:

* \(\omega_{cp}\);
* \(N_{cp}\);
* archivos CSV con los resultados;
* gráficas de la distribución.

---

## `03_calculo_eDcp.ipynb`

Calcula las distribuciones energéticas

$$
eD_{cp}(\epsilon)
$$

y

$$
eD_{cp}(\epsilon').
$$

A partir de sus máximos se determinan:

$$
\omega_\epsilon,
\qquad
\omega_{\epsilon'},
\qquad
\Delta_{Dcp}.
$$

---

## `04_resultados_analisis.ipynb`

Reúne los resultados obtenidos para los diferentes materiales y permite realizar la comparación final entre

$$
\frac{1}{2}\Delta_{Dcp}
$$

y la brecha superconductora conocida

$$
\Delta(0).
$$

---

# 📈 Resultados

El proyecto genera diferentes tipos de gráficas para cada material:

1. Función de Eliashberg \(\alpha^2F(\omega)\).
2. Densidad de estados fonónica.
3. Densidad de estados electrónica.
4. Distribución de pares de Cooper \(D_{cp}(\omega)\).
5. Distribución \(eD_{cp}(\epsilon)\).
6. Distribución \(eD_{cp}(\epsilon')\).

Las figuras se almacenan en:

```text
results/figures/
```

Las tablas numéricas se almacenan en:

```text
results/tables/
```

---

# 🧪 Materiales estudiados

Actualmente se consideran cinco materiales superconductores:

| Material | Símbolo | \(T_c\) aproximada |
| -------- | ------: | -----------------: |
| Aluminio |      Al |             1.18 K |
| Mercurio |      Hg |             4.15 K |
| Niobio   |      Nb |             9.25 K |
| Tantalio |      Ta |        4.47–4.48 K |
| Plomo    |      Pb |             7.19 K |

Los valores de \(T_c\) utilizados en los cálculos deben mantenerse consistentes con los valores definidos en los notebooks.

---

# ⚠️ Consideraciones numéricas

La calidad de los resultados depende directamente de la resolución de los datos de entrada.

En particular, la resolución de la densidad de estados electrónica puede ser crítica cuando las estructuras relevantes se encuentran en una escala energética mucho menor que el espaciamiento entre los puntos originales.

### Caso particular: Pb

Para Pb se encontró que el espaciamiento de los datos electrónicos originales es aproximadamente

$$
50\ \mathrm{meV},
$$

mientras que la frecuencia de corte obtenida para el cálculo es aproximadamente

$$
\omega_c\approx9.87\ \mathrm{meV}.
$$

Por lo tanto, dentro de la ventana física

$$
[-\omega_c,\omega_c]
$$

no existe una cantidad suficiente de puntos electrónicos originales para resolver directamente la estructura energética superconductora.

En consecuencia, los resultados asociados a Pb deben interpretarse con precaución y no deben considerarse equivalentes en calidad numérica a materiales cuyos datos electrónicos poseen mayor resolución.

---

# 📌 Resultados principales

Los cálculos permiten obtener las siguientes magnitudes:

| Material | \(\omega_\epsilon\) (meV) | \(\omega_{\epsilon'}\) (meV) | \(\Delta_{Dcp}\) (meV) | \(\frac12\Delta_{Dcp}\) (meV) |
| -------- | ------------------------: | ---------------------------: | ---------------------: | ----------------------------: |
| Al       |                    0.1655 |                      -0.1655 |                 0.3309 |                        0.1655 |
| Hg       |                    0.4257 |                      -0.4257 |                 0.8514 |                        0.4257 |
| Nb       |                    0.9414 |                      -0.9414 |                 1.8827 |                        0.9414 |
| Ta       |                    0.5336 |                      -0.5336 |                 1.0673 |                        0.5336 |
| Pb       |                    1.0455 |                      -1.0620 |                 2.1075 |                        1.0537 |

Estos valores se utilizan posteriormente para comparar la escala energética obtenida mediante la distribución de pares de Cooper con la brecha superconductora \(\Delta(0)\).

> **Nota:** los valores numéricos dependen de los datos de entrada, el preprocesamiento y la resolución de las DOS. En particular, el resultado de Pb presenta una limitación de resolución que debe considerarse en la interpretación.

---

# 🔎 Interpretación

El cálculo permite estudiar si existe una correspondencia cuantitativa entre la energía característica de la distribución de pares de Cooper y la brecha superconductora.

La comparación se realiza mediante

$$
\frac12\Delta_{Dcp}
\quad\text{vs.}\quad
\Delta(0).
$$

Los resultados muestran que la correspondencia depende del material y de la calidad de los datos utilizados.

Por ejemplo, para Al:

$$
\frac12\Delta_{Dcp}
\approx
0.1655\ \mathrm{meV},
$$

mientras que

$$
\Delta(0)
\approx
0.170\ \mathrm{meV}.
$$

La proximidad entre ambos valores es notable.

Sin embargo, esta comparación **no debe interpretarse como una identidad universal**

$$
\boxed{
\Delta(0)=\frac12\Delta_{Dcp}
}
$$

para todos los superconductores.

La relación debe evaluarse material por material y teniendo en cuenta las aproximaciones numéricas y la resolución de los datos de entrada.

---

# 🛠️ Tecnologías utilizadas

El proyecto utiliza principalmente:

* **Python**
* **NumPy** — cálculo numérico y manipulación de arreglos.
* **SciPy** — integración e interpolación numérica.
* **Pandas** — procesamiento de datos tabulares.
* **Matplotlib** — visualización científica.
* **Jupyter Notebook** — desarrollo y análisis interactivo.

---

# 📚 Organización del código

### `src/io.py`

Funciones encargadas de cargar los diferentes formatos de datos y convertirlos a la representación utilizada internamente.

### `src/units.py`

Centraliza las conversiones de unidades.

### `src/preprocessing.py`

Contiene las funciones de limpieza, ordenamiento, recorte e interpolación de los datos.

### `src/plotting.py`

Contiene las funciones utilizadas para generar las figuras del análisis.

### `src/cooper_pairs.py`

Contiene las funciones relacionadas directamente con el cálculo numérico de las distribuciones de pares de Cooper.

---

# 📖 Referencias

El desarrollo teórico y computacional del proyecto se basa principalmente en la guía de trabajo:

> G. I. González-Pedreros, *Cooper pair distribution function*, guía de trabajo, 4 de septiembre de 2026.

Las referencias adicionales utilizadas para \(T_c\), \(\Delta(0)\), datos de Eliashberg y propiedades de los materiales deben registrarse en el informe final del proyecto.

---

# 👨‍🔬 Proyecto académico

Proyecto desarrollado como parte de un trabajo de **Física Computacional / Ciencia de Materiales**, orientado al estudio numérico de superconductores convencionales y la distribución energética de pares de Cooper.

**Autor:** Sergio David Monroy Barragán

**Área:** Física computacional — Materia condensada — Superconductividad

---

<p align="center">
  <i>From electronic and phononic density of states to Cooper-pair distributions.</i>
</p>
