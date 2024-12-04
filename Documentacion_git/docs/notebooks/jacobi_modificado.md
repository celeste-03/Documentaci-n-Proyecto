# Método de Relajación de Jacobi Modificado

Este notebook implementa una mejora al método de Jacobi, añadiendo sobre-relajación para acelerar la convergencia en sistemas lineales.

```python
import numpy as np
import matplotlib.pyplot as plt

def jacobi_relaxation_mod(N, tolerance, L, V1, V2, omega):
    # Inicializar la matriz phi (potencial eléctrico)
    phi = np.zeros((N + 1, N + 1), dtype=float)

    # Aplicar condiciones de frontera
    phi[int((N / L) * 2):int((N / L) * 8), int((N / L) * 2)] = V1
    phi[int((N / L) * 2):int((N / L) * 8), int((N / L) * 8)] = V2

    delta = 1.0
    iteration = 0

    while delta > tolerance:
        phi_new = phi.copy()
        iteration += 1
        phi_new[1:N, 1:N] = (1.0 + omega) * 0.25 * (
            phi[2:N + 1, 1:N] + phi[0:N - 1, 1:N] + 
            phi[1:N, 2:N + 1] + phi[1:N, 0:N - 1]
        ) - omega * phi[1:N, 1:N]

        phi_new[int((N / L) * 2):int((N / L) * 8), int((N / L) * 2)] = V1
        phi_new[int((N / L) * 2):int((N / L) * 8), int((N / L) * 8)] = V2

        delta = np.max(np.abs(phi - phi_new))
        phi = phi_new

    print(f"El método convergió en {iteration} iteraciones.")
    return phi, iteration

phi, num_iterations = jacobi_relaxation_mod(100, 1e-3, 10, 1, -1, 1.5)


El método convergió en 523 iteraciones.

/tmp/ipykernel_5055/757925980.py:23: RuntimeWarning: overflow encountered in add
  phi[2:N + 1, 1:N] + phi[0:N - 1, 1:N] +
/tmp/ipykernel_5055/757925980.py:25: RuntimeWarning: overflow encountered in multiply
  ) - omega * phi[1:N, 1:N]
/tmp/ipykernel_5055/757925980.py:22: RuntimeWarning: overflow encountered in subtract
  phi_new[1:N, 1:N] = (1.0 + omega) * 0.25 * (
/tmp/ipykernel_5055/757925980.py:32: RuntimeWarning: overflow encountered in subtract
  delta = np.max(np.abs(phi - phi_new))
/tmp/ipykernel_5055/757925980.py:23: RuntimeWarning: invalid value encountered in add
  phi[2:N + 1, 1:N] + phi[0:N - 1, 1:N] +

# Graficar el resultado después de la convergencia, en este caso diverge este método.
#plt.figure(figsize=(8, 6))
#plt.imshow(phi, extent=[0, 10, 0, 10], origin='lower', cmap='viridis', interpolation='bilinear')
#plt.colorbar(label='Potencial eléctrico (phi)')
#plt.title(f'Distribución del Potencial Eléctrico - Iteraciones: {num_iterations}')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.show()

```
