/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobenius_norm.h"
#include "utils.h"
#include "d_alloc3d.h"
#include <omp.h>

void jacobi_inner_off_2(double ***u, double ***old_u, double ***f, const int N) {

    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);
    
    #pragma omp target teams loop is_device_ptr(u, old_u, f) 
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                
                u[i][j][k] = one_sixth * (old_u[i-1][j][k] + old_u[i+1][j][k] +old_u[i][j-1][k] + old_u[i][j+1][k] + old_u[i][j][k-1] + old_u[i][j][k+1] + grid_spacing_sq * f[i][j][k]);
            }
        }
    }
}

int jacobi_off_2(double *** u, double *** old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius) {
    
    
    
    int iter = 0;

    double* dev_u_data;
    double*** dev_u = d_malloc_3d(N+2, N+2, N+2, &dev_u_data);
    
    double* dev_u_old_data;
    double*** dev_u_old = d_malloc_3d(N+2, N+2, N+2, &dev_u_old_data);

    double* dev_f_data;
    double*** dev_f = d_malloc_3d(N+2, N+2, N+2, &dev_f_data);


    omp_target_memcpy(u[0][0], dev_u_data, sizeof(double)*(N+2)*(N+2)*(N+2), 0, 0, omp_get_default_device(), device(0));
    
    omp_target_memcpy(f[0][0],dev_f_data, sizeof(double)*(N+2)*(N+2)*(N+2), 0, 0, omp_get_default_device(), omp_get_initial_device());
    
    
    double delta_norm = INFINITY;
    
       while (delta_norm > threshold && iter < iter_max) {
        double ***tmp = dev_u;
        dev_u = dev_u_old;
        dev_u_old = tmp;
        jacobi_inner_off_2(dev_u, dev_u_old, dev_f, N);
        /*if (frobenius) {
            delta_norm = frobenius_norm(u, old_u, N);
        }*/
        ++iter;
    }

    omp_target_memcpy(dev_u_data, u, (N+2)*(N+2)*(N+2), 0, 0, omp_get_initial_device(), omp_get_default_device());

    d_free_3d(dev_u, dev_u_data);
    d_free_3d(dev_u_old, dev_u_old_data);
    d_free_3d(dev_f, dev_f_data);

    return iter;

}