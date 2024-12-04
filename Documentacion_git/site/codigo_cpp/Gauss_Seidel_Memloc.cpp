#include <iostream>
#include <vector>
#include <cmath> // Para std::abs
#include <algorithm> // Para std::max
#include <omp.h> //Open MP
#include <sys/time.h>//Para seconds()

double seconds()
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  
  return sec;
}

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

    // Se convierte las dimensiones físicas del capacitor a índices en la grilla
    int start_y = (N/L) * 2.0;  // Coordenada inicial (y = 2cm)
    int end_y = (N/L) * 8.0;    // Coordenada final (y= 8cm)
    int plate_x1 = (N / L) * 2.0; // Línea vertical izquierda (x = 2cm)
    int plate_x2 = (N / L) * 8.0; // Línea vertical derecha (x = 8cm)

    // Condiciones iniciales en las placas
    for (int i = start_y; i < end_y; ++i) {
        phi[i * grilla_size + plate_x1] = 1.0;  // +1V en la placa izquierda
        phi[i * grilla_size + plate_x2] = -1.0; // -1V en la placa derecha
    }

    std::vector<double> phi_copy = phi; // Copia inicial de phi

    double delta = 1.0;
    int its = 0;
    int num_procs;
    double time_1 = seconds();
   // omp_set_num_threads(3);
    while (delta > tolerance) {
        its++;
#pragma omp parallel
	{
	num_procs = omp_get_num_threads();
        // Se actualiza los valores de phi en toda la grilla excepto en los bordes
#pragma omp for
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
		    if((i>=start_y && i<=end_y)&&(j==plate_x1||j==plate_x2)){
            }
		    else{
                           phi[i * grilla_size + j] = (1.0 / 4.0) *
                           (phi[(i + 1) * grilla_size + j] +
                           phi[(i - 1) * grilla_size + j] +
                           phi[i * grilla_size + (j + 1)] +
                           phi[i * grilla_size + (j - 1)]);

		    }
	    }
        }

        // Restauro el potencial en las placas ya que este sí fue modificado
        /*for (int i = start_y; i < end_y; ++i) {
            phi[i * grilla_size + plate_x1] = 1.0;
            phi[i * grilla_size + plate_x2] = -1.0;
        }*/
	}

        delta = calcular_delta(phi, phi_copy);

        phi_copy = phi; // Actualizamos phi_copy con los nuevos valores de phi
	
}
    double time_2 = seconds();
    std::cout << "Se llego a la tolerancia tras " << its << " iteraciones." << std::endl;
    std::cout << "Number of Threads: " << num_procs << std::endl;
    std::cout << "Time to complete loop: " << time_2 - time_1 << std::endl;
    return 0;
}
