/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>
#include "utils.h"

double
gauss_seidel_inner_ser(double ***u, double ***f, const int N) {
    double total_delta = 0;

    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);

    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){

                const double old_u_value = u[i][j][k];

                u[i][j][k] = one_sixth * (u[i-1][j][k] + u[i+1][j][k] + u[i][j-1][k] + u[i][j+1][k] + u[i][j][k-1] + u[i][j][k+1] + grid_spacing_sq * f[i][j][k]);

                const double delta = old_u_value - u[i][j][k];
                total_delta += delta * delta;
            }
        }
    }
    return sqrt(total_delta);
}

int
gauss_seidel_ser(double ***u, double ***f, const int N, const int iter_max, const double threshold) {
    int iter = 0;
    double delta_norm = INFINITY;

    while (delta_norm > threshold && iter < iter_max) {
        delta_norm = gauss_seidel_inner_ser(u, f, N);
        ++iter;
    }

    return iter;
}
