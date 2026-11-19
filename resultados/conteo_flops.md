# Conteo de operaciones

Convención: 1 flop = 1 operación aritmética elemental (+, −, ×, ÷) en doble precisión.
Las comparaciones, los valores absolutos y los intercambios de fila no son flops: se
reportan aparte como **sobrecosto de estabilidad**. Los contadores se incrementan dentro
de los bucles de `src/gauss.py`; ninguna cifra de este archivo está escrita a mano.

Incluye: actualización de b durante la eliminación y sustitución regresiva.
No incluye: el refinamiento iterativo ni la verificación (se cuentan en `main.py`).

### Sin pivoteo

| Paso | Mult. | Div. | Sumas/restas | Total flops | Comparaciones | Div. de escala | Intercambios |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| k = 1 | 132 | 11 | 132 | 275 | 0 | 0 | 0 |
| k = 2 | 110 | 10 | 110 | 230 | 0 | 0 | 0 |
| k = 3 | 90 | 9 | 90 | 189 | 0 | 0 | 0 |
| k = 4 | 72 | 8 | 72 | 152 | 0 | 0 | 0 |
| k = 5 | 56 | 7 | 56 | 119 | 0 | 0 | 0 |
| k = 6 | 42 | 6 | 42 | 90 | 0 | 0 | 0 |
| k = 7 | 30 | 5 | 30 | 65 | 0 | 0 | 0 |
| k = 8 | 20 | 4 | 20 | 44 | 0 | 0 | 0 |
| k = 9 | 12 | 3 | 12 | 27 | 0 | 0 | 0 |
| k = 10 | 6 | 2 | 6 | 14 | 0 | 0 | 0 |
| k = 11 | 2 | 1 | 2 | 5 | 0 | 0 | 0 |
| sust. regresiva | 66 | 12 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **638** | **78** | **638** | **1354** | **0** | **0** | **0** |

### Pivoteo parcial

| Paso | Mult. | Div. | Sumas/restas | Total flops | Comparaciones | Div. de escala | Intercambios |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| k = 1 | 132 | 11 | 132 | 275 | 11 | 0 | 1 |
| k = 2 | 110 | 10 | 110 | 230 | 10 | 0 | 1 |
| k = 3 | 90 | 9 | 90 | 189 | 9 | 0 | 1 |
| k = 4 | 72 | 8 | 72 | 152 | 8 | 0 | 1 |
| k = 5 | 56 | 7 | 56 | 119 | 7 | 0 | 1 |
| k = 6 | 42 | 6 | 42 | 90 | 6 | 0 | 1 |
| k = 7 | 30 | 5 | 30 | 65 | 5 | 0 | 0 |
| k = 8 | 20 | 4 | 20 | 44 | 4 | 0 | 1 |
| k = 9 | 12 | 3 | 12 | 27 | 3 | 0 | 1 |
| k = 10 | 6 | 2 | 6 | 14 | 2 | 0 | 0 |
| k = 11 | 2 | 1 | 2 | 5 | 1 | 0 | 1 |
| sust. regresiva | 66 | 12 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **638** | **78** | **638** | **1354** | **66** | **0** | **9** |

### Pivoteo parcial escalado

| Paso | Mult. | Div. | Sumas/restas | Total flops | Comparaciones | Div. de escala | Intercambios |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| escalas de fila | 0 | 0 | 0 | 0 | 132 | 0 | 0 |
| k = 1 | 132 | 11 | 132 | 275 | 11 | 12 | 1 |
| k = 2 | 110 | 10 | 110 | 230 | 10 | 11 | 0 |
| k = 3 | 90 | 9 | 90 | 189 | 9 | 10 | 1 |
| k = 4 | 72 | 8 | 72 | 152 | 8 | 9 | 1 |
| k = 5 | 56 | 7 | 56 | 119 | 7 | 8 | 1 |
| k = 6 | 42 | 6 | 42 | 90 | 6 | 7 | 1 |
| k = 7 | 30 | 5 | 30 | 65 | 5 | 6 | 1 |
| k = 8 | 20 | 4 | 20 | 44 | 4 | 5 | 1 |
| k = 9 | 12 | 3 | 12 | 27 | 3 | 4 | 1 |
| k = 10 | 6 | 2 | 6 | 14 | 2 | 3 | 0 |
| k = 11 | 2 | 1 | 2 | 5 | 1 | 2 | 0 |
| sust. regresiva | 66 | 12 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **638** | **78** | **638** | **1354** | **198** | **77** | **8** |

## Comparación con la cota de referencia

| Estrategia | Flops medidos | Cota del enunciado | Desviación |
| --- | ---: | ---: | ---: |
| Sin pivoteo | 1354 | ≈ 1486 | -8.9 % |
| Pivoteo parcial | 1354 | ≈ 1486 | -8.9 % |
| Pivoteo parcial escalado | 1354 | ≈ 1486 | -8.9 % |

La cota ≈ 1486 del enunciado equivale a calcular también la columna k en cada paso
(a(i,k) − m·a(k,k), valor que se sabe nulo): son n(n−1) = 132 flops de más. El código
no los ejecuta —la subdiagonal se fija en cero por construcción y el multiplicador queda
guardado en su lugar—, por lo que el conteo medido (1354) queda 8,9 % bajo la cota,
dentro del ±10 % aceptado.

## Costo del refinamiento (aparte, no entra en el conteo de T1)

Cada corrección d = A⁻¹ r se resuelve con los factores L y U ya calculados (sin volver a eliminar): **276 flops por pasada** en doble precisión
(progresiva + regresiva). El residuo se evalúa en aritmética racional exacta y no se
cuenta como flop de binary64.

Resumen de la corrida escalada: factorización + sust. regresiva = 1354 flops; sobrecosto de estabilidad = 198 comparaciones + 77 divisiones de escala + 8 intercambios de fila.
