#!/usr/bin/env python


import numpy as np

def gauss_seidel(N=100, tolerance=1e-5, L=10, V1=1, V2=-1):
    # Iniciación de phi
    phi = np.zeros((N + 1, N + 1), dtype=float)
    
    # Condiciones iniciales (placas) (establecer el borde en 0V sería redundante)
    phi[int((N/L)*2):int((N/L)*8),int((N/L)*2)] = V1
    phi[int((N/L)*2):int((N/L)*8),int((N/L)*8)] = V2
    
    # Copia para evaluar el error
    phi_copy = phi.copy()
    
    delta = 1.0
    its = 0

    while delta > tolerance:
        its += 1
        
        # Itero. Aca debo usar for loops, para aseguarme de que los accesos se realicen entrada por entrada
        # de forma ordenada, ya que Gauss-Seidel depende de la actualización inmediata de los valores intermedios,
        # y el array slicing implica una actualizacion simultanea (util pero en Jacobi).
        for i in range(1,N):
            for j in range(1,N):
                phi[i, j] = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])/4.0
        

        # Restauro el potencial en las placas ya que lo modifiqué
        phi[int((N / L) * 2):int((N / L) * 8), int((N / L) * 2)] = V1
        phi[int((N / L) * 2):int((N / L) * 8), int((N / L) * 8)] = V2

        # Con delta se evalua si la solución se estabiliza
        delta = np.max(np.abs(phi - phi_copy))
        
        # Copiamos los valores de la nueva iteración
        phi_copy = phi.copy()

    return phi, its

gauss_vals, iterations = gauss_seidel()
print(iterations)


'''
En el jupyter:
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.imshow(gauss_vals)
plt.gray()
plt.show()
'''
