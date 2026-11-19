"""Eliminación gaussiana con y sin pivoteo, con contador de operaciones.

Todo el cálculo numérico se hace en doble precisión IEEE 754 (float de Python),
con bucles explícitos para que cada operación aritmética pase por el contador.

Modos de pivoteo:
    "sin"       -> sin intercambios (pivote = a(k,k) siempre)
    "parcial"   -> pivoteo parcial clásico (máximo |a(i,k)|)
    "escalado"  -> pivoteo parcial escalado por filas (máximo |a(i,k)| / s_i,
                   con s_i = máximo |a(i,j)| de la fila i en la matriz original)

La factorización es única: PA = LU. Los multiplicadores quedan guardados bajo
la diagonal de LU y se reutilizan para resolver otros lados derechos
(refinamiento iterativo) sin volver a eliminar.
"""

import math

from conteo import Contador


def _bits(x):
    """log2 de un valor positivo (instrumentación; no cuenta como flop)."""
    return math.log2(x)


def factorizar(A, b, modo="escalado", contador=None, registrar_traza=True):
    """Triangulariza [A | b] y resuelve por sustitución regresiva.

    Devuelve un diccionario con:
        x, LU (multiplicadores bajo la diagonal), perm, pasos (métricas por paso),
        trazas (matriz tras cada paso), crecimiento, mult_max, contador.
    """
    if modo not in ("sin", "parcial", "escalado"):
        raise ValueError("modo debe ser 'sin', 'parcial' o 'escalado'")
    n = len(A)
    c = contador if contador is not None else Contador()
    LU = [fila[:] for fila in A]
    y = list(b)
    perm = list(range(n))
    max_original = max(abs(v) for fila in A for v in fila)

    # Escalas de fila: se calculan UNA vez sobre la matriz original.
    escalas = None
    if modo == "escalado":
        c.paso("escalas")
        escalas = []
        for i in range(n):
            s = abs(LU[i][0])
            for j in range(1, n):
                c.cmp()
                if abs(LU[i][j]) > s:
                    s = abs(LU[i][j])
            escalas.append(s)

    pasos = []
    trazas = []
    crecimiento = 1.0
    mult_max_global = 0.0
    bits_abs_global = 0.0
    bits_canc_global = 0.0

    for k in range(n - 1):
        c.paso(k + 1)
        info = {"k": k + 1, "pivote_fila_original": None, "intercambio": False}

        # ---- búsqueda de pivote (sobrecosto de estabilidad) ----
        if modo != "sin":
            mejor = k
            if modo == "escalado":
                c.div_pivote()
                mejor_val = abs(LU[k][k]) / escalas[k]
            else:
                mejor_val = abs(LU[k][k])
            for i in range(k + 1, n):
                if modo == "escalado":
                    c.div_pivote()
                    val = abs(LU[i][k]) / escalas[i]
                else:
                    val = abs(LU[i][k])
                c.cmp()
                if val > mejor_val:
                    mejor_val = val
                    mejor = i
            if mejor != k:
                LU[k], LU[mejor] = LU[mejor], LU[k]
                y[k], y[mejor] = y[mejor], y[k]
                perm[k], perm[mejor] = perm[mejor], perm[k]
                if escalas is not None:
                    escalas[k], escalas[mejor] = escalas[mejor], escalas[k]
                c.intercambio()
                info["intercambio"] = True
        info["pivote_fila_original"] = perm[k] + 1
        piv = LU[k][k]
        info["pivote"] = piv
        if piv == 0.0:
            raise ZeroDivisionError(f"pivote nulo en el paso {k + 1}")

        # ---- eliminación: un multiplicador por fila, columnas k+1..n-1 ----
        mult_max = 0.0
        bits_abs = 0.0
        bits_canc = 0.0
        for i in range(k + 1, n):
            m = LU[i][k] / piv          # una sola división por fila
            c.div()
            LU[i][k] = m                # se guarda el multiplicador (L)
            if abs(m) > mult_max:
                mult_max = abs(m)
            for j in range(k + 1, n):   # la columna k no se toca: queda 0 por construcción
                t = m * LU[k][j]
                c.mul()
                a = LU[i][j]
                r = a - t
                c.add()
                LU[i][j] = r
                # instrumentación de mantisa (no es flop)
                if t != 0.0 and a != 0.0:
                    bits_abs = max(bits_abs, max(0.0, _bits(abs(t)) - _bits(abs(a))))
                if r != 0.0 and (a != 0.0 or t != 0.0):
                    mayor = max(abs(a), abs(t))
                    bits_canc = max(bits_canc, max(0.0, _bits(mayor) - _bits(abs(r))))
            t = m * y[k]
            c.mul()
            y[i] = y[i] - t
            c.add()

        # ---- métricas del paso ----
        # máximo |elemento| de la matriz de trabajo (U ya fija + submatriz activa)
        maximo = 0.0
        for i in range(n):
            j0 = i if i <= k else k + 1
            for j in range(j0, n):
                maximo = max(maximo, abs(LU[i][j]))
        crecimiento = max(crecimiento, maximo / max_original)
        mult_max_global = max(mult_max_global, mult_max)
        bits_abs_global = max(bits_abs_global, bits_abs)
        bits_canc_global = max(bits_canc_global, bits_canc)
        info.update(
            mult_max=mult_max,
            max_elemento=maximo,
            bits_absorbidos=bits_abs,
            bits_cancelados=bits_canc,
        )
        pasos.append(info)
        if registrar_traza:
            traza = []
            for i in range(n):
                fila = []
                for j in range(n):
                    fila.append(0.0 if (j <= k and i > j) else LU[i][j])
                traza.append(fila + [y[i]])
            trazas.append(traza)

    # ---- sustitución regresiva ----
    c.paso("regresiva")
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for j in range(i + 1, n):
            t = LU[i][j] * x[j]
            c.mul()
            s = s - t
            c.add()
        if LU[i][i] == 0.0:
            raise ZeroDivisionError("pivote nulo en la sustitución regresiva")
        x[i] = s / LU[i][i]
        c.div()

    return {
        "modo": modo,
        "x": x,
        "LU": LU,
        "perm": perm,
        "pasos": pasos,
        "trazas": trazas,
        "crecimiento": crecimiento,
        "mult_max": mult_max_global,
        "bits_absorbidos": bits_abs_global,
        "bits_cancelados": bits_canc_global,
        "contador": c,
        "y": y,
    }


def resolver_con_factores(LU, perm, r, contador=None):
    """Resuelve A d = r reutilizando PA = LU (sin volver a eliminar)."""
    n = len(LU)
    c = contador
    y = [r[perm[i]] for i in range(n)]
    for i in range(1, n):             # sustitución progresiva con L (diagonal unitaria)
        s = y[i]
        for j in range(i):
            t = LU[i][j] * y[j]
            s = s - t
            if c:
                c.mul()
                c.add()
        y[i] = s
    d = [0.0] * n
    for i in range(n - 1, -1, -1):    # sustitución regresiva con U
        s = y[i]
        for j in range(i + 1, n):
            t = LU[i][j] * d[j]
            s = s - t
            if c:
                c.mul()
                c.add()
        d[i] = s / LU[i][i]
        if c:
            c.div()
    return d


def resolver_transpuesta_con_factores(LU, perm, v):
    """Resuelve A^T w = v con PA = LU (A^T = U^T L^T P). Sirve para estimar κ."""
    n = len(LU)
    z = list(v)
    for i in range(n):                # U^T z = v  (U^T es triangular inferior)
        s = z[i]
        for j in range(i):
            s -= LU[j][i] * z[j]
        z[i] = s / LU[i][i]
    for i in range(n - 1, -1, -1):    # L^T u = z  (L^T es triangular superior, diag. 1)
        s = z[i]
        for j in range(i + 1, n):
            s -= LU[j][i] * z[j]
        z[i] = s
    w = [0.0] * n                     # w = P^T u
    for i in range(n):
        w[perm[i]] = z[i]
    return w
