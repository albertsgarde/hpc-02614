/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobenius_norm.h"
#include "utils.h"

void jacobi_inner_gpu_map(double ***u, double ***old_u, double ***f, const int N) {

    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);

    #pragma omp target teams loop map(from: u[0:N+2][0:N+2][0:N+2]) map(to: old_u[0:N+2][0:N+2][0:N+2]) map(to: f[0:N+2][0:N+2][0:N+2]) \
        collapse(2) num_teams(108) thread_limit(128)
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                double value1 = grid_spacing_sq * f[i][j][k];
                double value = one_sixth * (old_u[i-1][j][k] + old_u[i+1][j][k] + old_u[i][j-1][k] + old_u[i][j+1][k] 
                        + old_u[i][j][k-1] + old_u[i][j][k+1] + value1);
                u[i][j][k] = value;
            }
        }
    }
}

int jacobi_gpu_map(double *** u, double *** old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius) {
    int iter = 0;
    double delta_norm = INFINITY;
    
    /* #pragma omp target enter data \
        map(to: u[0:N][0:N][0:N]) map(alloc: old_u[0:N][0:N][0:N]) map(to: f[0:N][0:N][0:N]) */

    while (delta_norm > threshold && iter < iter_max) {
        double ***tmp = u;
        u = old_u;
        old_u = tmp;
        jacobi_inner_gpu_map(u, old_u, f, N);
        if (frobenius) {
            delta_norm = frobenius_norm(u, old_u, N);
        }
        ++iter;
    }
    
    /* #pragma omp target update from(u[0:N][0:N][0:N])
    #pragma omp target exit data \
        map(release: u[0:N][0:N][0:N]) map(release: old_u[0:N][0:N][0:N]) map(release: f[0:N][0:N][0:N]) */

    return iter;
}
