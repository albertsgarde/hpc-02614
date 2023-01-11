/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobenius_norm.h"
#include "utils.h"

double 
jacobi_inner(double ***u, double ***old_u, double ***f, const int N) {
    double total_delta = 0;

    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);

    #pragma omp for
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                u[i][j][k] = one_sixth * (old_u[i-1][j][k] + old_u[i+1][j][k] +old_u[i][j-1][k] + old_u[i][j+1][k] + old_u[i][j][k-1] + old_u[i][j][k+1] + grid_spacing_sq * f[i][j][k]);

                const double delta = u[i][j][k] - old_u[i][j][k];
                total_delta += delta * delta;
            }
        }
    }
    return total_delta;
}

int jacobi(double *** u, double *** old_u, double ***f, const int N, const int iter_max, const double threshold) {
    int iter = 0;
    double delta_norm = INFINITY;
    double delta = 0.0;

    #pragma omp parallel shared(iter, delta_norm, delta, u, old_u, f, N, iter_max, threshold)
    {
    while (delta_norm > threshold && iter < iter_max) {
        #pragma omp master
        {
        double ***tmp = u;
        u = old_u; // Add note that this caused problems only when thread num was 2 or sometimes on 3.
        old_u = tmp;
        }
        #pragma omp barrier
        double _delta = jacobi_inner(u, old_u, f, N);

        #pragma omp atomic
        delta += _delta;

        #pragma omp master
        {
        delta_norm = sqrt(delta);
        //printf("delta_norm = %f\n", delta_norm);
        delta = 0.;
        ++iter;
        }
        #pragma omp barrier
    }
    }
    return iter;
}
