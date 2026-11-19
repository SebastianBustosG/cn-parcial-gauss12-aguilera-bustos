# El pivote: eliminación gaussiana con mínimo costo y máxima exactitud

## 1. Identificación

| Campo | Detalle |
| --- | --- |
| Estudiantes | Cristobal Aguilera - Sebastián Bustos |
| Sección | 2 |
| Asignatura | Computación Numérica (INF-321) · Ingeniería Civil Informática, UCM |
| Docentes | Juan Manuel Cabezas Contardo · Pedro Pesante |
| Fecha | Talca, 19 de noviembre de 2026 |
| Repositorio | https://github.com/SebastianBustosG/cn-parcial-gauss12-aguilera-bustos.git |

**Declaración de uso de asistentes de IA.** Se usó Claude (Anthropic) como asistente para: leer el enunciado, proponer la arquitectura del repositorio, escribir una primera versión del código de `src/` y del borrador de este informe, y revisar la coherencia entre los números reportados y las salidas de `main.py`. Todas las cifras de este informe provienen de la ejecución del código entregado. [El/la estudiante declara que revisó el código línea por línea y que puede explicar la estrategia de eliminación, el conteo y el análisis de error en la defensa oral.]

## 2. Planteamiento

Se resuelve A x = b con A de 12×12 densa (144 elementos no nulos, 9 cifras significativas por elemento) y b de 18 a 20 cifras significativas. Ambos se leyeron sin truncar desde `datos/` y se validó que b conserva todas sus cifras. A se convierte a binary64 con redondeo correcto; el residuo y la solución de referencia usan los decimales exactos (`fractions.Fraction`).

Lo que hace costoso el caso:

- **Escalas de fila incompatibles.** El máximo |a(i,j)| por fila va de 9.3·10⁻⁵ (F1) a 9.7·10⁵ (F6): diez órdenes de magnitud. Restar un múltiplo de una fila grande a una pequeña desplaza los bits bajos del operando pequeño fuera de los 53 bits de mantisa.
- **Pivote inicial minúsculo.** a(1,1) = 1.3407·10⁻⁸ está doce órdenes bajo a(2,1) = 1.19·10⁴.
- **Casi dependencia lineal.** La fila 7 es casi combinación lineal de las filas 3 y 5 (sección 6). κ₂(A) ≈ 3.08·10¹³.

Ningún pivote se anula, así que el sistema se resuelve incluso sin pivoteo; por eso las corridas son comparables.

## 3. Estrategia de eliminación

Se factoriza **una sola vez** PA = LU sobre la matriz aumentada [A | b]. En cada paso k = 1…11:

1. (Con pivoteo) se elige la fila pivote y se intercambia. Las comparaciones y las divisiones de escala se cuentan aparte.
2. Para cada fila i > k se calcula **un solo multiplicador** m(i,k) = a(i,k)/a(k,k) y se reutiliza en toda la fila, incluido el lado derecho.
3. El barrido va de la columna k+1 en adelante: la columna k no se recalcula, porque se sabe nula. En su lugar se guarda m(i,k), de modo que L queda almacenada bajo la diagonal.
4. Al terminar, sustitución regresiva sobre U.

No se calcula la inversa de A. El refinamiento iterativo y la estimación de κ reutilizan L y U, sin factorizar de nuevo. Esto evita las dos fuentes de exceso que menciona el enunciado: recalcular multiplicadores por columna (≈ 1750 flops) y pasar por la inversa (≈ 7000 o más).

Se ejecutaron tres corridas del mismo código (`gauss.factorizar`): sin pivoteo, pivoteo parcial y pivoteo parcial escalado. En el escalado, s(i) = máx |a(i,j)| se calcula una vez sobre la matriz original y se elige la fila de mayor |a(i,k)|/s(i).

## 4. Conteo de operaciones

**Convención.** 1 flop = una suma, resta, multiplicación o división en binary64. Las comparaciones, los valores absolutos y los intercambios de fila no son flops: se reportan como sobrecosto de estabilidad. Las divisiones |a(i,k)|/s(i) del pivoteo escalado también se reportan aparte. Los contadores se incrementan dentro de los bucles de `src/gauss.py`; las tablas completas de las tres corridas están en `resultados/conteo_flops.md`.

Desglose por paso (idéntico para las tres corridas; incluye la actualización de b):

| Paso | Mult. | Div. | Sumas/restas | Total |
| --- | ---: | ---: | ---: | ---: |
| k = 1 | 132 | 11 | 132 | 275 |
| k = 2 | 110 | 10 | 110 | 230 |
| k = 3 | 90 | 9 | 90 | 189 |
| k = 4 | 72 | 8 | 72 | 152 |
| k = 5 | 56 | 7 | 56 | 119 |
| k = 6 | 42 | 6 | 42 | 90 |
| k = 7 | 30 | 5 | 30 | 65 |
| k = 8 | 20 | 4 | 20 | 44 |
| k = 9 | 12 | 3 | 12 | 27 |
| k = 10 | 6 | 2 | 6 | 14 |
| k = 11 | 2 | 1 | 2 | 5 |
| Sust. regresiva | 66 | 12 | 66 | 144 |
| **Total** | **638** | **78** | **638** | **1354** |

Multiplicaciones/divisiones: 716. Sumas/restas: 638.

**Comparación con la cota.** El término cúbico es 2n³/3 = 1152; la eliminación de A sola cuesta 1078 (= 2n³/3 − n²/2 − n/6), a lo que se suman 132 flops de actualizar b y 144 de la sustitución regresiva. La cota ≈ 1486 del enunciado corresponde a calcular además la columna k en cada paso (132 flops sobre valores que se sabe nulos). El código no los ejecuta, por lo que el total medido, 1354, queda 8.9 % bajo la cota: dentro del ±10 % aceptado y sin trabajo redundante.

**Sobrecosto del pivoteo (aparte):**

| Corrida | Comparaciones | Divisiones de escala | Intercambios | Flops base |
| --- | ---: | ---: | ---: | ---: |
| Sin pivoteo | 0 | 0 | 0 | 1354 |
| Parcial | 66 | 0 | 9 | 1354 |
| Parcial escalado | 198 (132 de escalas + 66) | 77 | 8 | 1354 |

El pivoteo no agrega ninguna operación aritmética a la eliminación; solo agrega comparaciones y, en la variante escalada, 77 divisiones. Cada pasada de refinamiento cuesta 276 flops en doble precisión (progresiva y regresiva con L y U).

## 5. Análisis de mantisa y pivoteo

Resultados de las tres corridas (`resultados/traza_eliminacion.txt`):

| Corrida | Mult. máx | Factor de crecimiento | Bits absorbidos | Bits cancelados | Error rel. ‖·‖∞ | Cifras correctas (peor componente) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sin pivoteo | 2.75·10¹³ | 2.64·10³ | 18.1 | 14.3 | 2.08·10⁻¹⁰ | 7.0 |
| Parcial | 6.2·10⁻¹ | 1.00 | 7.3 | 13.5 | 5.52·10⁻¹³ | 9.6 |
| Parcial escalado | 2.47·10⁹ | 1.78 | 6.3 | 6.3 | 2.36·10⁻¹³ | 10.0 |

Definiciones: el factor de crecimiento es el máximo |a(i,j)| de la matriz de trabajo en todos los pasos dividido por el máximo original (9.67·10⁵). Los bits absorbidos son el máximo de log₂(|m·a(k,j)|/|a(i,j)|), es decir, cuántos bits bajos de a(i,j) quedan por debajo del último bit de la suma. Los bits cancelados son el máximo de log₂(máx(|a|,|m·a|)/|a − m·a|). Ambas son métricas de instrumentación propias, aproximadas; el error relativo se mide contra la solución exacta del sistema.

**Por qué a(1,1) es inaceptable.** Con pivote 1.3407·10⁻⁸ el multiplicador de la fila 6 es 3.69·10⁵/1.3407·10⁻⁸ = 2.75·10¹³ ≈ 2⁴⁴·⁶. Tras el primer paso las filas 2 a 12 quedan dominadas por −m·a(1,j): el máximo de la matriz sube de 9.7·10⁵ a 2.55·10⁹ (factor 2.64·10³) y, en las filas de escala pequeña, el sumando original a(i,j) queda hasta 2¹⁸·¹ ≈ 2.8·10⁵ veces por debajo del término restado. De los 53 bits de esa suma, 18.1 dejan de representar a(i,j): sobreviven ≈ 35 bits (≈ 10.5 cifras) en el peor elemento, y los pasos siguientes (cancelaciones de hasta 14 bits) los reducen más. El orden de magnitud es compatible con lo observado: 7 cifras correctas al final. Sin pivoteo se pierden 3 cifras que con pivoteo se conservan.

**Pivoteo parcial frente a escalado.** El parcial clásico mantiene el crecimiento en 1.0 y multiplicadores menores que 1, pero elige pivotes por tamaño absoluto, sin mirar la escala de cada fila; el escalado compara cada candidato contra su propia fila. Aquí el escalado reduce a la mitad los bits absorbidos y cancelados (6.3 frente a 7.3 y 13.5), pero la diferencia final es pequeña (10.0 frente a 9.6 cifras), porque una vez elegidos pivotes razonables el error queda gobernado por el condicionamiento de la matriz con filas equilibradas (sección 6). Notar también que el escalado produce un multiplicador máximo de 2.47·10⁹ en el paso 3 sin dañar nada: el pivote (8.8·10⁻⁵, fila F1) es pequeño en valor absoluto pero grande respecto de su propia fila, y el crecimiento total sigue en 1.78. La medida relevante es el factor de crecimiento, no el multiplicador.

## 6. Condicionamiento

Estimaciones propias (sin invertir A, usando los factores L y U):

| Cantidad | Valor |
| --- | ---: |
| κ₁(A), estimador de Hager/Higham | 4.47·10¹³ |
| κ₂(A), potencias + iteración inversa | 3.085·10¹³ |
| κ₂(A), control `numpy.linalg.cond` | 3.085·10¹³ |
| σ_máx, σ_mín | 2.26·10⁶, 7.32·10⁻⁸ |
| κ₂·ε con ε = 2.22·10⁻¹⁶ | 6.85·10⁻³ |
| κ₂(A) con filas equilibradas | 8.38·10³ |
| κ₂(A_eq)·ε | 1.86·10⁻¹² |

**Cifras recuperables según la cota.** ‖δx‖/‖x‖ ≲ κ·ε ≈ 6.9·10⁻³ garantiza apenas 2.2 cifras. En la práctica se obtienen 10 (peor componente) y hasta 13.

**Por qué la cota es pesimista aquí.** κ₂(A) está inflado casi por completo por la dispersión de escalas de las filas, no por mal condicionamiento intrínseco. Al dividir cada fila por una potencia de 2 cercana a su máximo (escalado exacto en binary64), κ₂ cae a 8.4·10³ y la cota κ·ε pasa a 1.9·10⁻¹², frente al error observado de 2.4·10⁻¹³. El pivoteo escalado y el residuo por filas trabajan, en efecto, sobre el problema equilibrado. Además, la cota supone el peor alineamiento de los errores con la dirección de σ_mín.

**Casi dependencia de la fila 7.** Ajustando por mínimos cuadrados la fila 7 normalizada sobre las filas 3 y 5 normalizadas se obtiene F7 ≈ 0.989·F3 + 1.140·F5 con error relativo de ajuste 2.3·10⁻³ (el enunciado la describe como una perturbación de orden 10⁻⁴; la diferencia puede deberse a cómo se mide la perturbación y no cambia la conclusión). La confirmación está en el vector singular izquierdo de σ_mín de la matriz equilibrada: pesos F7 = 0.79, F5 = 0.45, F3 = 0.41. En la matriz sin equilibrar ese vector apunta solo a F5 (peso 1.000), simplemente la fila de menor norma: el artefacto de escala tapa la dependencia real. Esa es la razón de equilibrar antes de diagnosticar.

## 7. Resultados

Solución de la eliminación con pivoteo parcial escalado, redondeada a las cifras afirmadas, y solución tras refinar (`resultados/solucion.csv`):

| i | x (escalado, sin refinar) | Cifras afirmadas | Cifras medidas | x refinada (17 cifras) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.738491562700 | 12 | 12 | 0.73849156270000000 |
| 2 | −1.20593741680 | 12 | 12 | −1.2059374168000000 |
| 3 | 0.004387195200 | 10 | 10 | 0.0043871952000000000 |
| 4 | −0.091726453800 | 11 | 11 | −0.091726453800000000 |
| 5 | 2.846109372500 | 13 | 13 | 2.8461093725000000 |
| 6 | −0.00628154090 | 9 | 9 | −0.0062815409000000000 |
| 7 | 1.473820591600 | 13 | 13 | 1.4738205916000000 |
| 8 | −0.038492716400 | 11 | 11 | −0.038492716400000000 |
| 9 | 0.915264803700 | 12 | 12 | 0.91526480370000000 |
| 10 | −2.607139485200 | 13 | 13 | −2.6071394852000000 |
| 11 | 0.051873629400 | 11 | 11 | 0.051873629400000000 |
| 12 | −1.084937162500 | 13 | 13 | −1.0849371625000000 |

**Cómo se defienden las cifras.** Para la solución sin refinar, el número de cifras afirmadas sale de la primera corrección del refinamiento: |d(i)| estima el error de x(i), y se afirman p = ⌊E + 1 − log₁₀(2|d(i)|)⌋ cifras, con E el exponente decimal de x(i) (error menor que media unidad de la p-ésima cifra). Esta estimación no usa la solución verdadera; luego se contrastó contra la solución exacta del sistema (calculada con Fractions) y en las 12 incógnitas coinciden las cifras afirmadas con las medidas. Las componentes pequeñas (x₃, x₆) conservan menos cifras relativas: 10 y 9. Tras dos pasadas de refinamiento con residuo exacto, el error baja a ~10⁻³⁸ y se reportan 17 cifras, el máximo con sentido para este formato; los valores coinciden con la solución exacta del sistema entregado.

## 8. Verificación

**Residuo en precisión extendida.** Con r = b − A x calculado en `Fraction` (A y b exactos), la norma infinito para la solución escalada sin refinar es **2.29·10⁻¹⁰** (máximo en la fila 6, la de mayor escala). El mismo residuo en doble precisión da 5.8·10⁻¹⁰, porque suma términos de hasta ~10⁵ y queda afectado por su propio redondeo.

**Refinamiento iterativo** (residuo exacto, corrección con los factores L y U ya calculados):

| Pasada | ‖r‖∞ | Error rel. ‖·‖∞ | Cifras (peor componente) |
| ---: | ---: | ---: | ---: |
| 0 (sin refinar) | 2.29·10⁻¹⁰ | 2.36·10⁻¹³ | 10.0 |
| 1 | 1.37·10⁻²² | 6.2·10⁻²⁶ | 22.6 |
| 2 | 3.93·10⁻³⁵ | 1.6·10⁻³⁸ | 35.1 |

Cada pasada reduce el error en un factor ≈ 4·10¹², lejos de lo que sugeriría κ·ε ≈ 7·10⁻³. Como contraste, partiendo de la corrida sin pivoteo el refinamiento también converge (7.0 → 16.9 → 26.8 cifras en dos pasadas), pero el punto de partida es peor.

**Control con biblioteca.** `numpy.linalg.solve` (LAPACK) sobre la misma matriz en binary64 entrega entre 9.6 y 12.8 cifras correctas por componente, comparable a la implementación propia (10.0 a 13.2). Se usa solo como control.

## 9. Conclusiones

**El compromiso elegido.** Pivoteo parcial escalado por filas, una sola factorización, un multiplicador por fila, y refinamiento con residuo exacto. El costo aritmético de la eliminación es el mismo que sin pivoteo (1354 flops); el pivoteo agrega 198 comparaciones, 77 divisiones y 8 intercambios, y cada pasada de refinamiento agrega 276 flops. A cambio: factor de crecimiento 2.64·10³ → 1.78, y exactitud de 7 → 10 cifras antes de refinar y de 10 → 17 tras dos pasadas. No hay una tensión real entre costo y estabilidad en la eliminación: la estabilidad se compra con comparaciones y unas pocas divisiones, no con aritmética del núcleo.

**Qué haría distinto con n = 1000.** El costo base sería ≈ 2n³/3 ≈ 6.7·10⁸ flops, y el sobrecosto del pivoteo (O(n²) comparaciones) se vuelve despreciable. Usaría una LU por bloques (BLAS nivel 3) por eficiencia de caché; equilibraría las filas antes de factorizar; estimaría κ con Hager (O(n²) por solve) en lugar de cualquier cálculo de inversa; y haría el refinamiento con residuo en precisión extendida (suma compensada o doble-doble en vez de `Fraction`, que no escala), cada pasada O(n²) frente al O(n³) de la factorización. Si la memoria apremia, factorizar en precisión simple y refinar en doble es la opción clásica.

## 10. Reproducibilidad

- Python 3.12.3 (se requiere ≥ 3.10), numpy 2.5.4; el resto de la biblioteca estándar (`fractions`, `decimal`, `math`).
- Comando: `pip install -r requirements.txt` y luego `python src/main.py`. Regenera los cuatro archivos de `resultados/` en menos de un segundo.
- El cálculo es determinista: no hay números aleatorios.
- Fuentes de método consultadas: enunciado de la prueba (Computación Numérica, INF-321, UCM); N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2.ª ed., SIAM, 2002 (estimador de condición y refinamiento iterativo); G. H. Golub y C. F. Van Loan, *Matrix Computations*, 4.ª ed., Johns Hopkins, 2013 (eliminación con pivoteo y análisis de error).
