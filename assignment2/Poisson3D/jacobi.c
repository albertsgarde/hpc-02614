/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobeniusnorm.c"

void jacobi(double ***U, double ***oldU, double ***f, int N, int kmax, double threshold) {
    
    int k = 0;
    int d = 100000000;

    while (d > threshold && k < kmax) {

        oldU = U;
        jacobiInner(U, oldU, f, N);
        d = frobeniusnorm(U, oldU, (N+2));
        k += 1;
    }

}

void jacobiInner(double ***U, double ***oldU, double ***f, int N) {

    double one_sixth = 1/6;
    double gridsizeSq = 1/(N*N) ;
    
    for (int i = 1; i < (N + 2); i++){
        for (int j=1; j < (N + 2); j++){
            for (int k=1; k < (N + 2); k++){

                newval = (oldU[i-1][j][k] + oldU[i+1][j][k] +oldU[i][j-1][k] + oldU[i][j]+1[k] + oldU[i][j][k-1] + oldU[i][j][k+1])
  
                U[i][j][k] = one_sixth * (newval + gridsizeSq * f[i][j][k]);
            }
        }
    }


}
