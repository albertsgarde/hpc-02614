/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobenius_norm.c"

void jacobi_inner(double *const *const *const U, double *const *const *const old_U, const double *const *const *const f, const int N) {

    const double one_sixth = 1./6.;
    // See comments in gauss_seidel.c
    const double grid_size_sq = 1./(double)(N*N) ;
    
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                U[i][j][k] = one_sixth * (old_U[i-1][j][k] + old_U[i+1][j][k] +old_U[i][j-1][k] + old_U[i][j+1][k] + old_U[i][j][k-1] + old_U[i][j][k+1] + grid_size_sq * f[i][j][k]);
            }
        }
    }
}

void jacobi(double *const *const * U, double *const *const * old_U, const double *const *const *const f, const int N, const int iter_max, const double threshold) {
    int iter = 0;
    double delta_norm = INFINITY;

    while (delta_norm > threshold && iter < iter_max) {
        double *const *const *const tmp = U;
        U = old_U;
        old_U = tmp;
        jacobi_inner(U, old_U, f, N);
        delta_norm = frobenius_norm((const double *const *const *)U, (const double *const *const *)old_U, N);
        ++iter;
    }
}
