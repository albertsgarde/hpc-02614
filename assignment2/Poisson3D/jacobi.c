/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobenius_norm.h"

void jacobi_inner(double ***u, double ***old_u, double ***f, const int N) {

    const double one_sixth = 1./6.;
    // See comments in gauss_seidel.c
    const double grid_size_sq = 1./(double)(N*N) ;
    
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                u[i][j][k] = one_sixth * (old_u[i-1][j][k] + old_u[i+1][j][k] +old_u[i][j-1][k] + old_u[i][j+1][k] + old_u[i][j][k-1] + old_u[i][j][k+1] + grid_size_sq * f[i][j][k]);
            }
        }
    }
}

void jacobi(double *** u, double *** old_u, double ***f, const int N, const int iter_max, const double threshold) {
    int iter = 0;
    double delta_norm = INFINITY;

    while (delta_norm > threshold && iter < iter_max) {
        double ***tmp = u;
        u = old_u;
        old_u = tmp;
        jacobi_inner(u, old_u, f, N);
        delta_norm = frobenius_norm(u, old_u, N);
        ++iter;
    }
}
