
import math
import os
import platform
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import conteo  # noqa: E402
import gauss  # noqa: E402
import verificacion as ver  # noqa: E402
from conteo import Contador  # noqa: E402

RAIZ = os.path.dirname(AQUI)
DATOS = os.path.join(RAIZ, "datos")
SALIDA = os.path.join(RAIZ, "resultados")
getcontext().prec = 60
EPS = 2.0 ** -52          # épsilon de máquina de binary64 (≈ 2.22e-16)
DIGITOS_REFINADA = 17     # tope útil: binary64 y resolución de los datos (ver informe)


def fmt(v, cifras):
    """Formatea con exactamente `cifras` cifras significativas (nunca más)."""
    return f"{float(v):.{cifras - 1}e}"


def escribir(nombre, texto):
    with open(os.path.join(SALIDA, nombre), "w", encoding="utf-8") as f:
        f.write(texto)


def main():
    os.makedirs(SALIDA, exist_ok=True)

    # ------------------------------------------------------------------ datos
    A_txt, b_txt = ver.leer_datos(
        os.path.join(DATOS, "matriz_a.csv"), os.path.join(DATOS, "vector_b.csv")
    )
    ver.validar_datos(A_txt, b_txt)
    A_ex, b_ex = ver.a_fracciones(A_txt, b_txt)      # decimales exactos
    A, b = ver.a_floats(A_txt, b_txt)                # binary64 correctamente redondeado
    n = len(A)

    # Solución exacta del sistema tal como se entregó (referencia; Gauss propio con Fractions)
    x_exacta = ver.resolver_exacto(A_ex, b_ex)

    # ------------------------------------------------------------- 3 corridas
    corridas = {}
    for modo in ("sin", "parcial", "escalado"):
        corridas[modo] = gauss.factorizar(A, b, modo, Contador())

    def err_inf(x):
        num = max(abs(Fraction(xi) - r) for xi, r in zip(x, x_exacta))
        den = max(abs(r) for r in x_exacta)
        return float(num / den)

    def dig_min(x):
        return min(ver.digitos_correctos(xi, r) for xi, r in zip(x, x_exacta))

    nombres = {
        "sin": "Sin pivoteo",
        "parcial": "Pivoteo parcial",
        "escalado": "Pivoteo parcial escalado",
    }

    # ----------------------------------------------------------------- trazas
    t = ["TRAZA DE LA ELIMINACIÓN GAUSSIANA (n = 12)", "=" * 78, ""]
    t.append("RESUMEN COMPARATIVO DE LAS CORRIDAS (T2)")
    t.append("")
    t.append(
        f"{'corrida':<26}{'mult. máx':>12}{'factor crec.':>14}{'bits absorb.':>14}"
        f"{'bits cancel.':>14}{'err. rel. ∞':>14}{'díg. mín.':>11}"
    )
    for modo in ("sin", "parcial", "escalado"):
        r = corridas[modo]
        t.append(
            f"{nombres[modo]:<26}{r['mult_max']:>12.3e}{r['crecimiento']:>14.3e}"
            f"{r['bits_absorbidos']:>14.1f}{r['bits_cancelados']:>14.1f}"
            f"{err_inf(r['x']):>14.2e}{dig_min(r['x']):>11.1f}"
        )
    t += [
        "",
        "mult. máx   : máximo |m(i,k)| de toda la eliminación.",
        "factor crec.: máx |a(i,j)^(k)| sobre todos los pasos / máx |a(i,j)| original.",
        "bits absorb.: máx log2(|m·a(k,j)| / |a(i,j)|): bits bajos de a(i,j) que la suma deja",
        "              fuera de los 53 bits de mantisa (el sumando grande los desplaza).",
        "bits cancel.: máx log2(máx(|a|,|m·a|) / |a − m·a|): bits perdidos por cancelación.",
        "err. rel. ∞ : ||x − x_exacta||_inf / ||x_exacta||_inf (x_exacta: Fractions, ver informe).",
        "",
    ]
    for modo in ("sin", "escalado"):
        r = corridas[modo]
        t += ["-" * 78, f"CORRIDA: {nombres[modo].upper()}", "-" * 78, ""]
        for paso, traza in zip(r["pasos"], r["trazas"]):
            t.append(
                f"Paso k = {paso['k']:>2} | pivote = {paso['pivote']:+.9e} "
                f"(fila original F{paso['pivote_fila_original']}"
                f"{', intercambio' if paso['intercambio'] else ''}) | "
                f"mult. máx = {paso['mult_max']:.3e} | máx |elem| = {paso['max_elemento']:.3e}"
            )
            t.append("  matriz [A | b] tras el paso (ceros bajo la diagonal en columnas ≤ k):")
            for fila in traza:
                t.append("  " + " ".join(f"{v:>16.9e}" for v in fila))
            t.append("")
        bajo = sum(
            1
            for i in range(n)
            for j in range(i)
            if r["trazas"][-1][i][j] != 0.0
        )
        t.append(f"Elementos no nulos bajo la diagonal tras el paso 11: {bajo} (de 66)")
        t.append(f"Permutación de filas final (índices originales, base 1): "
                 f"{[p + 1 for p in r['perm']]}")
        t.append("")
    escribir("traza_eliminacion.txt", "\n".join(t))

    # --------------------------------------------------------- conteo de flops
    lista = [(nombres[m], corridas[m]) for m in ("sin", "parcial", "escalado")]
    texto_conteo = conteo.reporte(lista)

    # ------------------------------------------------- refinamiento iterativo
    esc = corridas["escalado"]
    c_ref = Contador()
    x_ref, hist = ver.refinar(A_ex, b_ex, esc["x"], esc["LU"], esc["perm"], 10, c_ref)
    c_ref_sin = Contador()
    _, hist_sin = ver.refinar(A_ex, b_ex, corridas["sin"]["x"], corridas["sin"]["LU"],
                              corridas["sin"]["perm"], 10, c_ref_sin)
    c1 = Contador()
    c1.paso("una pasada")
    gauss.resolver_con_factores(esc["LU"], esc["perm"], [1.0] * n, c1)
    flops_pasada = c1.flops_base()

    texto_conteo += "\n".join(
        [
            "",
            "## Costo del refinamiento (aparte, no entra en el conteo de T1)",
            "",
            f"Cada corrección d = A⁻¹ r se resuelve con los factores L y U ya calculados "
            f"(sin volver a eliminar): **{flops_pasada} flops por pasada** en doble precisión",
            "(progresiva + regresiva). El residuo se evalúa en aritmética racional exacta y no se",
            "cuenta como flop de binary64.",
            "",
            "Resumen de la corrida escalada: factorización + sust. regresiva = "
            f"{esc['contador'].flops_base()} flops; sobrecosto de estabilidad = "
            f"{esc['contador'].total('cmp')} comparaciones + "
            f"{esc['contador'].total('div_pivote')} divisiones de escala + "
            f"{esc['contador'].total('swap')} intercambios de fila.",
            "",
        ]
    )
    escribir("conteo_flops.md", texto_conteo)

    # ------------------------------------------------------ solución (T3)
    d0 = hist[0]["d"]
    claim_gauss = []
    for i in range(n):
        p_est = ver.digitos_estimados(esc["x"][i], d0[i])
        claim_gauss.append(max(1, min(17, math.floor(p_est))))
    filas_csv = [
        "i,x_gauss_escalado,digitos_afirmados_gauss,digitos_reales_gauss,"
        "x_refinada,digitos_afirmados_refinada,digitos_reales_refinada"
    ]
    avisos = []
    for i in range(n):
        txt_g = fmt(esc["x"][i], claim_gauss[i])
        d_ref = Decimal(x_ref[i].numerator) / Decimal(x_ref[i].denominator)
        txt_r = f"{d_ref:.{DIGITOS_REFINADA - 1}e}"
        # dígitos efectivamente correctos del valor calculado SIN redondear a lo afirmado
        real_g = ver.digitos_correctos(esc["x"][i], x_exacta[i])
        real_r = min(17.0, ver.digitos_correctos(x_ref[i], x_exacta[i]))
        if claim_gauss[i] > math.floor(real_g):
            avisos.append(f"x{i + 1}: se afirmaron {claim_gauss[i]} dígitos y hay {real_g:.1f}")
        if DIGITOS_REFINADA > math.floor(real_r):
            avisos.append(f"x{i + 1} refinada: se afirmaron {DIGITOS_REFINADA} y hay {real_r:.1f}")
        filas_csv.append(
            f"{i + 1},{txt_g},{claim_gauss[i]},{math.floor(real_g)},"
            f"{txt_r},{DIGITOS_REFINADA},{math.floor(real_r)}"
        )
    escribir("solucion.csv", "\n".join(filas_csv) + "\n")

    # --------------------------------------------- condicionamiento y control
    LUe, pe = esc["LU"], esc["perm"]
    resolver = lambda r: gauss.resolver_con_factores(LUe, pe, r)
    resolver_t = lambda v: gauss.resolver_transpuesta_con_factores(LUe, pe, v)
    cond = ver.condicionamiento(A, resolver, resolver_t)
    # Matriz con filas equilibradas (potencias de 2: el escalado es exacto en binary64).
    # Se reutilizan los mismos factores L y U: Ae = D^-1 A  =>  Ae^-1 y = A^-1 (D y).
    esc2 = [2.0 ** round(math.log2(max(abs(v) for v in fila))) for fila in A]
    Ae = [[v / si for v in fila] for fila, si in zip(A, esc2)]
    resolver_e = lambda y: resolver([d * si for d, si in zip(y, esc2)])
    resolver_te = lambda v: [si * w for si, w in zip(esc2, resolver_t(v))]
    cond_e = ver.condicionamiento(Ae, resolver_e, resolver_te)
    x_np, k2_np, det_np = ver.control_numpy(A, b)
    (c3, c5), err7 = ver.casi_dependencia_fila(A, 6, (2, 4))
    izq = cond["vec_izq_sigma_min"]
    orden_izq = sorted(range(n), key=lambda i: -abs(izq[i]))[:3]
    izq_e = cond_e["vec_izq_sigma_min"]
    orden_izq_e = sorted(range(n), key=lambda i: -abs(izq_e[i]))[:3]
    cota = cond["kappa2_est"] * EPS
    cota_e = cond_e["kappa2_est"] * EPS

    # ------------------------------------------------------------ residuo.txt
    r0 = ver.residuo_exacto(A_ex, b_ex, esc["x"])
    r_doble = [
        float(b[i]) - sum(A[i][j] * esc["x"][j] for j in range(n)) for i in range(n)
    ]
    L = ["RESIDUO, REFINAMIENTO Y CONDICIONAMIENTO", "=" * 78, ""]
    L += ["1. Residuo r = b − A·x en precisión extendida (Fractions; A y b exactos)", ""]
    L.append("x = solución de la eliminación con pivoteo parcial escalado (binary64)")
    for i in range(n):
        L.append(f"  r[{i + 1:>2}] = {float(r0[i]):+.6e}")
    L.append(f"  ||r||_inf (extendida) = {norma_inf_(r0):.6e}")
    L.append(f"  ||r||_inf (doble, para contraste) = {max(abs(v) for v in r_doble):.6e}")
    L.append("  (el residuo en doble precisión suma términos de hasta ~1e5 y queda afectado por su propio")
    L.append("   redondeo: solo el residuo en precisión extendida sirve para juzgar la exactitud)")
    L += ["", "2. Refinamiento iterativo (residuo exacto, corrección con los factores L y U)", ""]
    L.append(f"{'iter':>4} {'||r||_inf':>14} {'||d||/||x||':>14} {'err. rel. ∞ vs exacta':>24} {'díg. mín.':>10}")
    for h in hist:
        xi = h["x"]
        L.append(
            f"{h['iter']:>4} {h['norma_r']:>14.3e} {h['rel_d']:>14.3e} "
            f"{err_inf(xi):>24.3e} {dig_min(xi):>10.1f}"
        )
    L += ["", "   Contraste: el mismo refinamiento partiendo de la corrida SIN pivoteo", ""]
    L.append(f"{'iter':>4} {'||r||_inf':>14} {'||d||/||x||':>14} {'err. rel. ∞ vs exacta':>24} {'díg. mín.':>10}")
    for h in hist_sin:
        xi = h["x"]
        L.append(
            f"{h['iter']:>4} {h['norma_r']:>14.3e} {h['rel_d']:>14.3e} "
            f"{err_inf(xi):>24.3e} {dig_min(xi):>10.1f}"
        )
    L += ["", f"   Residuo final (escalado, tras refinar): ||r||_inf = "
              f"{norma_inf_(ver.residuo_exacto(A_ex, b_ex, x_ref)):.3e}"]
    L += ["", "3. Condicionamiento", ""]
    L.append(f"  κ1(A)  estimado (Hager/Higham sobre L y U) = {cond['kappa1_est']:.3e}")
    L.append(f"  κ2(A)  estimado (potencias + iteración inversa) = {cond['kappa2_est']:.3e}")
    L.append(f"        σ_max ≈ {cond['sigma_max']:.3e}   σ_min ≈ {cond['sigma_min']:.3e}")
    L.append(f"  Control numpy.linalg.cond(A, 2) = {k2_np:.3e}   det(A) = {det_np:.4e}")
    L.append(f"  Cota clásica ||δx||/||x|| ≲ κ2·ε = {cota:.2e}   (ε = {EPS:.3e})")
    L.append(f"  Dígitos que garantiza la cota ≈ {max(0.0, -math.log10(cota)):.1f}")
    L.append("")
    L.append("  Con las filas equilibradas (cada fila dividida por una potencia de 2 cercana a su máximo):")
    L.append(f"    κ2(A_eq) estimado = {cond_e['kappa2_est']:.3e}   κ2·ε = {cota_e:.2e}")
    L.append(f"    Dígitos que garantiza esa cota ≈ {max(0.0, -math.log10(cota_e)):.1f}")
    L.append("")
    L.append("  Casi dependencia lineal (fila 7 frente a filas 3 y 5):")
    L.append(f"    F7/||F7|| ≈ {c3:+.6f}·F3/||F3|| {c5:+.6f}·F5/||F5||   (mínimos cuadrados, filas normalizadas)")
    L.append(f"    error relativo del ajuste = {err7:.2e}")
    L.append("    Vector singular izquierdo de σ_min, matriz original   : "
             + ", ".join(f"F{i + 1} ({abs(izq[i]):.3f})" for i in orden_izq))
    L.append("    Vector singular izquierdo de σ_min, filas equilibradas: "
             + ", ".join(f"F{i + 1} ({abs(izq_e[i]):.3f})" for i in orden_izq_e))
    L += ["", "4. Control con biblioteca (numpy.linalg.solve; solo verificación)", ""]
    L.append(f"{'i':>3} {'dígitos correctos numpy':>26} {'dígitos propios (escalado)':>28}")
    for i in range(n):
        L.append(
            f"{i + 1:>3} {ver.digitos_correctos(x_np[i], x_exacta[i]):>26.1f} "
            f"{ver.digitos_correctos(esc['x'][i], x_exacta[i]):>28.1f}"
        )
    L.append("")
    escribir("residuo.txt", "\n".join(L))

    # -------------------------------------------------------------- consola
    print("== Reproducción completa: resultados/ actualizado ==")
    print(f"Python {platform.python_version()}")
    for modo in ("sin", "parcial", "escalado"):
        r = corridas[modo]
        print(f"{nombres[modo]:<26} flops={r['contador'].flops_base():>5} "
              f"mult.máx={r['mult_max']:.3e} crec={r['crecimiento']:.3e} "
              f"err={err_inf(r['x']):.2e} díg.mín={dig_min(r['x']):.1f}")
    print(f"Refinamiento: ||r||∞ {hist[0]['norma_r']:.2e} -> {hist[-1]['norma_r']:.2e} "
          f"en {hist[-1]['iter']} pasadas")
    print(f"κ2 estimado = {cond['kappa2_est']:.3e} (numpy: {k2_np:.3e})")
    if avisos:
        print("AVISOS de dígitos afirmados:", *avisos, sep="\n  ")
    else:
        print("Dígitos afirmados coherentes con la solución exacta de referencia.")


def norma_inf_(v):
    return max(abs(float(t)) for t in v)


if __name__ == "__main__":
    main()
