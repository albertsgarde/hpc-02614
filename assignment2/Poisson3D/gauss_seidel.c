/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void
gauss_seidel(double ***U, double ***f, int N, int kmax, double threshold) {
    
    int k = 0;
    int d = 100000000;

    while (d > threshold && k < kmax) {

        double old_u_val = 0;
        double old_sum = 0;

        double one_sixth = 1/6;
        double gridsizeSq = 1/(N*N) ;
    
        for (int i = 1; i < (N + 2); i++){
            for (int j=1; j < (N + 2); j++){
                for (int k=1; k < (N + 2); k++){

                    old_u_val = U[i][j][k];

                    newval = (U[i-1][j][k] + U[i+1][j][k] + U[i][j-1][k] + U[i][j]+1[k] + U[i][j][k-1] + U[i][j][k+1])
    
                    U[i][j][k] = one_sixth * (newval + gridsizeSq * f[i][j][k]);

                    old_sum += old_u_val * old_u_val;
                }
            }
        }


        d = Math.sqrt(old_sum);
        k += 1;
    }

    




}
