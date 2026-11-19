"""Verificación independiente: residuo exacto, refinamiento iterativo, κ(A) y control.

La precisión extendida se hace con `fractions.Fraction` (aritmética racional exacta):
A y b se leen como decimales exactos (sin pasar por float), de modo que el residuo
r = b − A·x se calcula sin ningún error de redondeo.
"""

import csv
import math
from fractions import Fraction

from gauss import resolver_con_factores, resolver_transpuesta_con_factores


# ----------------------------------------------------------------------------
# Lectura de datos exactos
# ----------------------------------------------------------------------------
def leer_datos(ruta_a, ruta_b):
    """Devuelve (A_texto, b_texto): cadenas tal cual vienen en los CSV (sin truncar)."""
    with open(ruta_a, newline="", encoding="utf-8") as f:
        A_txt = [fila for fila in csv.reader(f) if fila]
    with open(ruta_b, encoding="utf-8") as f:
        b_txt = [linea.strip() for linea in f if linea.strip()]
    return A_txt, b_txt


def validar_datos(A_txt, b_txt):
    assert len(A_txt) == 12 and all(len(f) == 12 for f in A_txt), "A debe ser 12x12"
    assert len(b_txt) == 12, "b debe tener 12 componentes"
    cifras_b = [len(t.lstrip("-").replace(".", "").lstrip("0")) for t in b_txt]
    assert min(cifras_b) >= 18, "b fue truncado: se esperan 18 a 20 cifras significativas"


def a_fracciones(A_txt, b_txt):
    A = [[Fraction(t) for t in fila] for fila in A_txt]
    b = [Fraction(t) for t in b_txt]
    return A, b


def a_floats(A_txt, b_txt):
    """Lectura correctamente redondeada a binary64 (lo mejor que puede hacer un float)."""
    return [[float(t) for t in fila] for fila in A_txt], [float(t) for t in b_txt]


# ----------------------------------------------------------------------------
# Solución exacta de referencia (eliminación propia sobre racionales)
# ----------------------------------------------------------------------------
def resolver_exacto(A, b):
    """Gauss con Fractions (propio, no es rutina de biblioteca). Solo como referencia."""
    n = len(A)
    M = [fila[:] + [b[i]] for i, fila in enumerate(A)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= m * M[k][j]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


# ----------------------------------------------------------------------------
# Residuo en precisión extendida
# ----------------------------------------------------------------------------
def residuo_exacto(A, b, x):
    """r = b − A x exacto. x puede ser float o Fraction."""
    xf = [v if isinstance(v, Fraction) else Fraction(v) for v in x]
    n = len(A)
    return [b[i] - sum(A[i][j] * xf[j] for j in range(n)) for i in range(n)]


def norma_inf(v):
    return max(abs(float(t)) for t in v)


def digitos_correctos(x, ref):
    """Cifras significativas correctas de x respecto de ref.

    Convención: x tiene p cifras correctas si |x − ref| ≤ 0.5·10^(E − p + 1),
    con E = floor(log10|ref|) (error menor que media unidad de la p-ésima cifra).
    Devuelve p como número real; el valor entero defendible es su parte entera.
    """
    if ref == 0:
        return 0.0
    err = abs(Fraction(x) - ref)
    if err == 0:
        return 99.0
    E = math.floor(math.log10(abs(float(ref))))
    return E + 1 - math.log10(2.0 * float(err))


def digitos_estimados(x, d):
    """Mismo conteo, pero usando la corrección d (estimador del error) en lugar de ref."""
    if d == 0:
        return 99.0
    E = math.floor(math.log10(abs(float(x))))
    return E + 1 - math.log10(2.0 * abs(float(d)))


# ----------------------------------------------------------------------------
# Refinamiento iterativo (residuo exacto, correcciones con los factores L y U)
# ----------------------------------------------------------------------------
def refinar(A_exacta, b_exacta, x0, LU, perm, max_iter=10, contador=None):
    """Refinamiento con residuo en precisión extendida.

    Devuelve (x_final como Fraction, historial). Cada fila del historial:
    (iter, ||r||_inf, ||d||_inf / ||x||_inf, d_i) donde d es la corrección calculada.
    """
    if contador is not None:
        contador.paso("refinamiento")
    x = [Fraction(v) for v in x0]
    historial = []
    for it in range(max_iter + 1):
        r = residuo_exacto(A_exacta, b_exacta, x)
        rf = [float(t) for t in r]
        d = resolver_con_factores(LU, perm, rf, contador)
        nx = max(abs(float(v)) for v in x)
        historial.append(
            {
                "iter": it,
                "norma_r": norma_inf(r),
                "rel_d": max(abs(t) for t in d) / nx,
                "d": d,
                "x": list(x),
            }
        )
        if it == max_iter or historial[-1]["rel_d"] < 1e-30:
            break
        x = [xi + Fraction(di) for xi, di in zip(x, d)]
    return x, historial


# ----------------------------------------------------------------------------
# Condicionamiento
# ----------------------------------------------------------------------------
def _norma1(A):
    n = len(A)
    return max(sum(abs(A[i][j]) for i in range(n)) for j in range(n))


def _norma_inf_mat(A):
    return max(sum(abs(v) for v in fila) for fila in A)


def estimar_norma1_inversa(resolver, resolver_t, n, iteraciones=5):
    """Estimador de Hager/Higham de ||A^-1||_1 usando solo solves (no se forma la inversa).
    `resolver(r)` resuelve A d = r y `resolver_t(v)` resuelve A^T w = v."""
    x = [1.0 / n] * n
    est = 0.0
    for _ in range(iteraciones):
        y = resolver(x)
        est_nuevo = sum(abs(t) for t in y)
        xi = [1.0 if t >= 0 else -1.0 for t in y]
        z = resolver_t(xi)
        j = max(range(n), key=lambda i: abs(z[i]))
        if est_nuevo <= est or abs(z[j]) <= sum(zi * xi_ for zi, xi_ in zip(z, x)):
            est = max(est, est_nuevo)
            break
        est = est_nuevo
        x = [0.0] * n
        x[j] = 1.0
    return est


def estimar_kappa2(A, resolver, resolver_t, iteraciones=60):
    """κ2 ≈ σ_max/σ_min por iteración de potencias (σ_max) e iteración inversa (σ_min),
    esta última con solves de A y A^T sobre los factores (sin invertir A).
    Devuelve (κ2, σ_max, σ_min, vector singular izquierdo asociado a σ_min)."""
    n = len(A)

    def matvec(v):
        return [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]

    def matvec_t(v):
        return [sum(A[i][j] * v[i] for i in range(n)) for j in range(n)]

    def normalizar(v):
        nv = math.sqrt(sum(t * t for t in v))
        return [t / nv for t in v], nv

    v, _ = normalizar([1.0 + 0.1 * i for i in range(n)])
    smax = 0.0
    for _ in range(iteraciones):
        v, nw = normalizar(matvec_t(matvec(v)))
        smax = math.sqrt(nw)
    u, _ = normalizar([1.0 + 0.1 * i for i in range(n)])
    inv_smin = 0.0
    for _ in range(iteraciones):
        # (A^T A)^-1 u = A^-1 (A^-T u)
        u, nw = normalizar(resolver(resolver_t(u)))
        inv_smin = math.sqrt(nw)
    # vector singular izquierdo de σ_min: dirección de A^-T aplicada a u
    izq, _ = normalizar(resolver_t(u))
    return smax * inv_smin, smax, 1.0 / inv_smin, izq


def condicionamiento(A, resolver, resolver_t):
    n = len(A)
    n1 = _norma1(A)
    inv1 = estimar_norma1_inversa(resolver, resolver_t, n)
    k2, smax, smin, izq = estimar_kappa2(A, resolver, resolver_t)
    return {
        "kappa1_est": n1 * inv1,
        "kappa2_est": k2,
        "sigma_max": smax,
        "sigma_min": smin,
        "vec_izq_sigma_min": izq,
    }


def casi_dependencia_fila(A, objetivo, otras):
    """Ajusta (en mínimos cuadrados, 2 variables, ecuaciones normales 2x2) la fila
    `objetivo` como combinación de las filas `otras`. Devuelve (coef, error relativo).
    Se normaliza cada fila por su norma para evitar el problema de escalas."""
    n = len(A)

    def nrm(f):
        return math.sqrt(sum(t * t for t in f))

    fo = [t / nrm(A[objetivo]) for t in A[objetivo]]
    f1 = [t / nrm(A[otras[0]]) for t in A[otras[0]]]
    f2 = [t / nrm(A[otras[1]]) for t in A[otras[1]]]
    g11 = sum(a * a for a in f1)
    g12 = sum(a * b for a, b in zip(f1, f2))
    g22 = sum(b * b for b in f2)
    h1 = sum(a * c for a, c in zip(f1, fo))
    h2 = sum(b * c for b, c in zip(f2, fo))
    det = g11 * g22 - g12 * g12
    c1 = (h1 * g22 - h2 * g12) / det
    c2 = (g11 * h2 - g12 * h1) / det
    res = [fo[j] - c1 * f1[j] - c2 * f2[j] for j in range(n)]
    return (c1, c2), nrm(res) / nrm(fo)


# ----------------------------------------------------------------------------
# Control con biblioteca (NO sustituye la implementación propia)
# ----------------------------------------------------------------------------
def control_numpy(A_float, b_float):
    import numpy as np

    Af = np.array(A_float)
    bf = np.array(b_float)
    x = np.linalg.solve(Af, bf)
    return [float(t) for t in x], float(np.linalg.cond(Af, 2)), float(np.linalg.det(Af))
