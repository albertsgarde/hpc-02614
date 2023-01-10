/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void
gauss_seidel(double ***U, double ***f, int N, int kmax, double threshold) {
    
    int k = 0;
    int d = 100000000;

    while (d > threshold && k < kmax) {

        U = U;
        gauss_seidel_inner(U, U, f, N);
        d = frobeniusnorm(U, U, (N+2));
        k += 1;
    }

}

void gauss_seidel_inner(double ***U, double ***U, double ***f, int N) {

    double one_sixth = 1/6;
    double gridsizeSq = 1/(N*N) ;
    
    for (int i = 1; i < (N + 2); i++){
        for (int j=1; j < (N + 2); j++){
            for (int k=1; k < (N + 2); k++){

                newval = (U[i-1][j][k] + U[i+1][j][k] + U[i][j-1][k] + U[i][j]+1[k] + U[i][j][k-1] + U[i][j][k+1])
  
                U[i][j][k] = one_sixth * (newval + gridsizeSq * f[i][j][k]);
            }
        }
    }

}
