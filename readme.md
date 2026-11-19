# cn-parcial-gauss12 · El pivote

Resuelve un sistema lineal denso de 12 incógnitas por eliminación gaussiana propia (con y sin pivoteo), contando cada operación y verificando el resultado en precisión extendida.

**Autoría:** [Nombre Apellido] · Sección [N] · Computación Numérica (INF-321) · Ingeniería Civil Informática, Universidad Católica del Maule · Prueba Parcial, noviembre de 2026.

## Requisitos

- Python 3.10 o superior (probado con 3.12).
- Instalación de dependencias (solo `numpy`, usado para el control y para `cond`):

```bash
pip install -r requirements.txt
```

## Ejecución

```bash
python src/main.py
```

Un único comando, sin intervención manual, regenera todo `resultados/`:

| Archivo | Contenido |
| --- | --- |
| `resultados/solucion.csv` | las 12 incógnitas con los dígitos que se afirman correctos |
| `resultados/traza_eliminacion.txt` | resumen comparativo y matriz tras cada uno de los 11 pasos (sin pivoteo y con pivoteo escalado) |
| `resultados/conteo_flops.md` | operaciones por paso, convención y sobrecosto del pivoteo |
| `resultados/residuo.txt` | residuo exacto, refinamiento iterativo, κ(A) y control con numpy |

## Estructura

```
datos/        A (12x12) y b tal como se entregaron, sin truncar
src/          gauss.py (eliminación), conteo.py (flops), verificacion.py (residuo, κ), main.py
resultados/   salidas generadas por main.py
informe.md    informe técnico y análisis de error
```

## Resultado resumido

Pivoteo parcial escalado + dos pasadas de refinamiento con residuo exacto. Valores a 10 cifras significativas:

| i | x(i) | i | x(i) |
| ---: | ---: | ---: | ---: |
| 1 | 0.7384915627 | 7 | 1.473820592 |
| 2 | -1.205937417 | 8 | -0.03849271640 |
| 3 | 0.004387195200 | 9 | 0.9152648037 |
| 4 | -0.09172645380 | 10 | -2.607139485 |
| 5 | 2.846109373 | 11 | 0.05187362940 |
| 6 | -0.006281540900 | 12 | -1.084937163 |

- Costo: **1354 flops** (eliminación + actualización de b + sustitución regresiva), más 198 comparaciones, 77 divisiones de escala y 8 intercambios de fila.
- Residuo exacto ‖r‖∞: **2.29·10⁻¹⁰** antes de refinar, **1.37·10⁻²²** tras una pasada y **3.93·10⁻³⁵** tras dos.

El análisis completo está en [`informe.md`](informe.md).

## Licencia

Trabajo académico de uso educativo para la asignatura INF-321. No está pensado para redistribución.
