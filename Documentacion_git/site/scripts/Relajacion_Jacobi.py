#!/usr/bin/env python

import numpy as np
'''
Argumentos: 
    N: el tamaño lineal de la grilla cuadrada (NxN) 
    Tolerancia
    L: el tamaño lineal del capacitor cuadrado (LxL)
    Primer y segundo voltaje en las placas
'''
def jacobi_relaxation(N=100, tolerance=1e-5, L=10, V1=1, V2=-1):
    # Iniciación de phi
    phi = np.zeros((N + 1, N + 1), dtype=float)
    # Condiciones iniciales (placas) (establecer el borde en 0V sería redundante)
    phi[int((N/L)*2):int((N/L)*8),int((N/L)*2)] = V1
    phi[int((N/L)*2):int((N/L)*8),int((N/L)*8)] = V2

    # Acá en el codigo del profe se usa phiprime y luego se usa temp, esto se debe a que x=y hace shallow copies
    # por eso es necesario intercambiar phi y phiprime, para no trabajar sobre el MISMO espacio en memoria, yo voy a resolver esto 
    # haciendo un phinew, con un nombre mas descriptivo y que va a servir para un mejor uso de la memoria, ya que tan solo ocupo
    # almacenar la siguiente iteracion y compararla con la vieja.
    
    delta = 1.0
    its = 0
    
    while delta > tolerance:
        its += 1
        # Creo una copia (DEEP COPY) para el estado futuro
        phi_new = phi.copy()

        phi_new[1:N, 1:N] = (1.0/4.0) * (phi[2:N + 1, 1:N] + phi[0:N - 1, 1:N] + phi[1:N, 2:N + 1] + phi[1:N, 0:N - 1])
        #phiprime[1:N, 1:N] usa array slicing para implementar una vectorizacion con numpy (haciendolo mucho más eficiente 
        #que los for loops anidados), de manera que recorro toda la grilla, excepto los bordes ([1] hasta [N]). Ahora, si 
        #quiero phi[i+1, j] entonces se sumo 1 en la primera entrada -> phi[2:N + 1, 1:N]  

        # Restauro el potencial en las placas ya que este sí lo modifiqué
        phi_new[int((N / L) * 2):int((N / L) * 8), int((N / L) * 2)] = V1
        phi_new[int((N / L) * 2):int((N / L) * 8), int((N / L) * 8)] = V2
        
        # Con delta se evalua si la solución se estabiliza
        delta = np.max(np.abs(phi - phi_new))
        
        # Se guarda el valor de phiprime como el valor "presente" para la siguiente iteracion (phi)
        phi = phi_new
        
    return phi, its

jacobi_vals, iterations = jacobi_relaxation()
print(iterations)

'''
En el jupyter:
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.imshow(jacobi_vals)
plt.gray()
plt.show()
'''
