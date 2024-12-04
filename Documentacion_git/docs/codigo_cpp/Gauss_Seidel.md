# Solución de la Ecuación de Laplace en C++

Aqui se presenta una traducción del código implementado en python. Esto con el objetivo de mejorar su velocidad, y posteriormente paralelizar.

## Código Fuente

```cpp
#include <iostream>
#include <vector>
#include <cmath> // Para std::abs
#include <algorithm> // Para std::max

double calcular_delta(const std::vector<double> &phi, const std::vector<double> &phi_copy) {
    double delta = 0.0;
    for (size_t i = 0; i < phi.size(); ++i) {
        delta = std::max(delta, std::abs(phi[i] - phi_copy[i]));
    }
    return delta;
}

int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << " --N [Tamaño lineal de la grilla cuadrada (NxN)] --t [Tolerancia] --L [Tamaño lineal del capacitor cuadrado (LxL)]" << std::endl;
        return 1;
    }

    std::cout << "Bienvenido." << std::endl;
    std::cout << "Se esta calculando la solucion de la ecuacion de Laplace bajo los parametros definidos..." << std::endl;

    int N = atoi(argv[2]);  // Tamaño de la grilla (NxN)
    double tolerance = atof(argv[4]);  // Tolerancia
    int L = atoi(argv[6]);  // Tamaño lineal del capacitor fisico
    int grilla_size = N + 1;  // Tamaño total de la grilla (incluye frontera)

    std::vector<double> phi(grilla_size * grilla_size, 0.0);

    int start_y = (N/L) * 2.0;  // Coordenada inicial (y = 2cm)
    int end_y = (N/L) * 8.0;    // Coordenada final (y= 8cm)
    int plate_x1 = (N / L) * 2.0; // Línea vertical izquierda (x = 2cm)
    int plate_x2 = (N / L) * 8.0; // Línea vertical derecha (x = 8cm)

    for (int i = start_y; i < end_y; ++i) {
        phi[i * grilla_size + plate_x1] = 1.0;  // +1V en la placa izquierda
        phi[i * grilla_size + plate_x2] = -1.0; // -1V en la placa derecha
    }

    std::vector<double> phi_copy = phi;

    double delta = 1.0;
    int its = 0;

    while (delta > tolerance) {
        its++;

        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                phi[i * grilla_size + j] = (1.0 / 4.0) *
                    (phi[(i + 1) * grilla_size + j] +
                     phi[(i - 1) * grilla_size + j] +
                     phi[i * grilla_size + (j + 1)] +
                     phi[i * grilla_size + (j - 1)]);
            }
        }

        for (int i = start_y; i < end_y; ++i) {
            phi[i * grilla_size + plate_x1] = 1.0;
            phi[i * grilla_size + plate_x2] = -1.0;
        }

        delta = calcular_delta(phi, phi_copy);

        phi_copy = phi;
    }

    std::cout << "Se llego a la tolerancia tras " << its << " iteraciones." << std::endl;
    return 0;
}
```

Se logra la convergencia tras 1196 iteraciones (las mismas que presenta su contraparte en python).

