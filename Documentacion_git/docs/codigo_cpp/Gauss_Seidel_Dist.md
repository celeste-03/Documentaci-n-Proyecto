# Solución de la Ecuación de Laplace en C++ en memoria compartida 

Este código pretende utilizar MPI para lograr una paralelización en memoria compartida. En resumen cada proceso se queda con una banda horizontal de la placa, esta banda esta compuesta por una región sobre la que itera, y los valores colindantes. Luego, los procesos con `rank` par se comunican los de `rank` impar para actualizar dichos valores colindantes e iterar de manera independiente.



## Código Fuente

```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>

double calcular_delta(const std::vector<double> &phi, const std::vector<double> &phi_copy) {
    double delta = 0.0;
    for (size_t i = 0; i < phi.size(); ++i) {
        delta = std::max(delta, std::abs(phi[i] - phi_copy[i]));
    }
    return delta;
}


void iteracion_local(std::vector<double> &phi, int local_rows, int grilla_size) {
    for (int i = 1; i <= local_rows; ++i) {
        for (int j = 1; j < grilla_size - 1; ++j) {
            phi[i * grilla_size + j] = 0.25 * (
                phi[(i + 1) * grilla_size + j] +  // Celda inferior (ya actualizada en esta iteración)
                phi[(i - 1) * grilla_size + j] +  // Celda superior (ya actualizada en esta iteración)
                phi[i * grilla_size + (j + 1)] +  
                phi[i * grilla_size + (j - 1)] );   
        }
    }
}


void establecer_placas(std::vector<double> &local_phi, int start_row, int end_row, int N, int L, int grilla_size) {
    int start_y = (N / L) * 2;
    int end_y = (N / L) * 8;
    int plate_x1 = (N / L) * 2;
    int plate_x2 = (N / L) * 8;
    // Verificamos si la región de filas que maneja este proceso intersecta con las coordenadas de las placas
    if (start_row <= end_y && end_row >= start_y) {
        // Actualizamos las filas correspondientes en la región de las placas
        for (int i = std::max(start_row, start_y); i < std::min(end_row, end_y); ++i) {
            local_phi[(i - start_row + 1) * grilla_size + plate_x1] = 1.0;  // +1V en la placa izquierda
            local_phi[(i - start_row + 1) * grilla_size + plate_x2] = -1.0; // -1V en la placa derecha
        }
    }
}


int main(int argc, char* argv[]) {
    
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 7) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " --N [Tamaño lineal de la grilla cuadrada (NxN)] --t [Tolerancia] --L [Tamaño lineal del capacitor cuadrado (LxL)]" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        std::cout << "Bienvenido." << std::endl;
        std::cout << "Se esta calculando la solucion de la ecuacion de Laplace bajo los parametros definidos..." << std::endl;
    }

    int N = atoi(argv[2]);
    double tolerance = atof(argv[4]);
    int L = atoi(argv[6]);
    int grilla_size = N + 1;

    double time_1 = MPI_Wtime();
    
    // Dividir la grilla en bandas horizontales
    int nlocal = N / size; // Filas por proceso //rows_per_proc 
    int rest = N % size; // Lo que sobra si el tamaño no es divisible entre los procesos disponibles //extra_rows 
    
    int local_rows = nlocal;  // Inicializamos con el valor básico
    if (rank < rest) {
        local_rows += 1;  // Si hay resto, los primeros procesos tiene una fila adicional
    }
    int start_row = rank * nlocal + std::min(rank, rest); // A los últimos procesos se les suma el resto, los primeros se van acomodando con el mismmo espaciado entre ellos
    int end_row = start_row + local_rows;
    int local_rows_con_vecinos = local_rows + 2; // Añadir las filas solapadas //paded_rows

    
    std::vector<double> local_phi(local_rows_con_vecinos * grilla_size, 0.0);
    establecer_placas(local_phi, start_row, end_row, N, L, grilla_size);

    std::vector<double> local_phi_copy = local_phi;

    double global_delta = 1.0;
    double local_delta = 1.0;
    int its = 0;

    while (global_delta > tolerance) {
        its++;

        // Fase 1: Procesos pares envían, procesos impares reciben
        if (rank % 2 == 0) { // Procesos pares envían
            if (rank < size - 1) { // Si hay un vecino abajo
                MPI_Send(&local_phi[local_rows * grilla_size], grilla_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
            if (rank > 0) { // Si hay un vecino arriba
                MPI_Send(&local_phi[grilla_size], grilla_size, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            }
        } else { // Procesos impares reciben
            if (rank > 0) { // Si hay un vecino arriba
                MPI_Recv(&local_phi[0], grilla_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank < size - 1) { // Si hay un vecino abajo
                MPI_Recv(&local_phi[(local_rows + 1) * grilla_size], grilla_size, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD); // Sincronización entre envíos y recepciones

        // Fase 2: Procesos impares iteran localmente
        if (rank % 2 == 1) {
            iteracion_local(local_phi, local_rows, grilla_size);
            establecer_placas(local_phi, start_row, end_row, N, L, grilla_size);
        }

        // Fase 3: Procesos impares envían, procesos pares reciben
        if (rank % 2 == 1) { // Procesos impares envían
            if (rank < size - 1) { // Si hay un vecino abajo
                MPI_Send(&local_phi[local_rows * grilla_size], grilla_size, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD);
            }
            if (rank > 0) { // Si hay un vecino arriba
                MPI_Send(&local_phi[grilla_size], grilla_size, MPI_DOUBLE, rank - 1, 3, MPI_COMM_WORLD);
            }
        } else { // Procesos pares reciben
            if (rank > 0) { // Si hay un vecino arriba
                MPI_Recv(&local_phi[0], grilla_size, MPI_DOUBLE, rank - 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank < size - 1) { // Si hay un vecino abajo
                MPI_Recv(&local_phi[(local_rows + 1) * grilla_size], grilla_size, MPI_DOUBLE, rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD); // Sincronización entre envíos y recepciones

        // Fase 4: Procesos pares actualizan e iteran
        if (rank % 2 == 0) {
            iteracion_local(local_phi, local_rows, grilla_size);
            establecer_placas(local_phi, start_row, end_row, N, L, grilla_size);
        }

        // Calcular delta local
        local_delta = calcular_delta(local_phi, local_phi_copy);

        // Reducir el delta máximo global
        MPI_Allreduce(&local_delta, &global_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        local_phi_copy = local_phi;
    }

    double time_2 = MPI_Wtime();
    
    if (rank == 0) {
        std::cout << "Se llegó a la tolerancia tras " << its << " iteraciones." << std::endl;
        std::cout << "Tiempo utilizado: " << time_2 - time_1 << std::endl;
    }

    MPI_Finalize();
    return 0;
}
```

Necesita 1199 iteraciones.

