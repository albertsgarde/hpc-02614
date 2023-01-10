/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

double
gauss_seidel_inner(double *const *const *const U, const double *const *const *const f, const int N) {
    double total_delta = 0;

    const double one_sixth = 1./6.;
    // Maybe replace N with N+2? or replace 1.^2 with 2.^2?
    const double grid_size_sq = 1./(double)(N*N);

    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){

                const double old_u_value = U[i][j][k];

                U[i][j][k] = one_sixth * (U[i-1][j][k] + U[i+1][j][k] + U[i][j-1][k] + U[i][j+1][k] + U[i][j][k-1] + U[i][j][k+1] + grid_size_sq * f[i][j][k]);

                const double delta = old_u_value - U[i][j][k];
                total_delta += delta * delta;
            }
        }
    }
    return sqrt(total_delta);
}

void
gauss_seidel(double *const *const *const U, const double *const *const *const f, const int N, const int iter_max, const double threshold) {
    
    int iter = 0;
    double delta_norm = INFINITY;

    while (delta_norm > threshold && iter < iter_max) {
        delta_norm = gauss_seidel_inner(U, f, N);
        ++iter;
    }
}


