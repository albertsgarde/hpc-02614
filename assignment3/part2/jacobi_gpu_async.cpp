/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include <stdio.h>
#include "frobenius_norm.h"
#include "utils.h"
#include "d_alloc3d.h"
#include <omp.h>
#include "warm_up.h"

void jacobi_inner_gpu_async(double ***d_u_0, double ***d_old_u_0, double ***d_f_0, double ***d_u_1, double ***d_old_u_1, double ***d_f_1, const int N) {
    
    const int n_half = (N+2)/2;
    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);

    #pragma omp target nowait teams loop is_device_ptr(d_u_0, d_old_u_0, d_f_0) device(0)\
        collapse(3)  thread_limit(128) 
    for (int i = 1; i < n_half-1; i++){
        for (int j=1; j < N +1 ; j++){
            for (int k=1; k < N + 1; k++){
                //printf("Loop1 iteration: (%d,%d,%d) \n", i, j, k);
                double value1 = grid_spacing_sq * d_f_0[i][j][k];
                double value = one_sixth * (d_old_u_0[i-1][j][k] + d_old_u_0[i+1][j][k] + d_old_u_0[i][j-1][k] + d_old_u_0[i][j+1][k] 
                        + d_old_u_0[i][j][k-1] + d_old_u_0[i][j][k+1] + value1);
                d_u_0[i][j][k] = value;
            }
        }
    }
    
    omp_set_default_device(0);
    #pragma omp target nowait teams loop is_device_ptr(d_u_0, d_old_u_0, d_f_0,  d_u_1, d_old_u_1, d_f_1)\
        collapse(2) thread_limit(128) 
    for (int j=1; j < (N + 1); j++){
        for (int k=1; k < (N + 1); k++){
            //printf("Mid: (%d,%d,%d) \n", n_half, j,k);
            double value1 = grid_spacing_sq * d_f_0[n_half-1][j][k];
    
            double value = one_sixth * (d_old_u_0[n_half-2][j][k] + d_old_u_1[0][j][k] + d_old_u_0[n_half-1][j-1][k] + d_old_u_0[n_half-1][j+1][k] 
                        + d_old_u_0[n_half-1][j][k-1] + d_old_u_0[n_half-1][j][k+1] + value1);
            d_u_0[n_half-1][j][k] = value;
        }
    }
    omp_set_default_device(1);
    #pragma omp target nowait teams loop is_device_ptr(d_u_0, d_old_u_0, d_f_0,  d_u_1, d_old_u_1, d_f_1)\
        collapse(2) thread_limit(128) 
    for (int j=1; j < (N + 1); j++){
        for (int k=1; k < (N + 1); k++){
            //printf("Mid: (%d,%d,%d) \n", n_half, j,k);
            double value1 = grid_spacing_sq * d_f_1[0][j][k];
    
            double value = one_sixth * (d_old_u_0[n_half-1][j][k] + d_old_u_1[1][j][k] + d_old_u_1[0][j-1][k] + d_old_u_1[0][j+1][k] 
                        + d_old_u_1[0][j][k-1] + d_old_u_1[0][j][k+1] + value1);
            d_u_1[0][j][k] = value;
        }
    }
    
    #pragma omp target nowait teams loop is_device_ptr(d_u_1, d_old_u_1, d_f_1) device(1) \
        collapse(3) thread_limit(128) 
    for (int i = 1; i < n_half; i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                double value1 = grid_spacing_sq * d_f_1[i][j][k];
                double value = one_sixth * (d_old_u_1[i-1][j][k] + d_old_u_1[i+1][j][k] + d_old_u_1[i][j-1][k] + d_old_u_1[i][j+1][k] 
                        + d_old_u_1[i][j][k-1] + d_old_u_1[i][j][k+1] + value1);
                d_u_1[i][j][k] = value;
                }
            }
        }
    
    #pragma omp taskwait

}

int jacobi_gpu_async(double *** u, double *** old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius) {
    

    cudaSetDevice(0);
    cudaDeviceEnablePeerAccess(1, 0); // (dev 1, future flag)
    cudaSetDevice(1);
    cudaDeviceEnablePeerAccess(0, 0);
    cudaSetDevice(0);

    const int grid_size = N+2;
    int iter = 0;
    double delta_norm = INFINITY;

    omp_set_default_device(0);
    warm_up();
    double* d_u_data_0;
    double*** d_u_0 = d_malloc_3d(grid_size, grid_size, grid_size, &d_u_data_0);
    double* d_old_u_data_0;
    double*** d_old_u_0 = d_malloc_3d(grid_size, grid_size, grid_size, &d_old_u_data_0);
    double* d_f_data_0;
    double*** d_f_0 = d_malloc_3d(grid_size, grid_size, grid_size, &d_f_data_0);

    omp_set_default_device(1);
    warm_up();
    double* d_u_data_1;
    double*** d_u_1 = d_malloc_3d(grid_size, grid_size, grid_size, &d_u_data_1);
    double* d_old_u_data_1;
    double*** d_old_u_1 = d_malloc_3d(grid_size, grid_size, grid_size, &d_old_u_data_1);
    double* d_f_data_1;
    double*** d_f_1 = d_malloc_3d(grid_size, grid_size, grid_size, &d_f_data_1);

    double* fdata[] = {d_f_data_0, d_f_data_1};
    double* udata[] = {d_u_data_0, d_u_data_1};
    
    #pragma omp parallel for
    for (int w = 0; w < 2; w++) {
        omp_set_default_device(w);
      
        omp_target_memcpy(fdata[w], f[w*(N/2)][0], (grid_size/2) * grid_size * grid_size * sizeof(double),
                0, 0, omp_get_default_device(), omp_get_initial_device());
    
        omp_target_memcpy(udata[w], u[w*(N/2)][0], (grid_size/2) * grid_size * grid_size * sizeof(double),
                0, 0, omp_get_default_device(), omp_get_initial_device());
    
    }


    while (delta_norm > threshold && iter < iter_max) {
   
        double ***tmp0 = d_old_u_0;
        d_old_u_0 = d_u_0;
        d_u_0 = tmp0;

        double ***tmp1 = d_old_u_1;
        d_old_u_1 = d_u_1;
        d_u_1 = tmp1;
        jacobi_inner_gpu_async(d_u_0, d_old_u_0, d_f_0, d_u_1, d_old_u_1, d_f_1, N);
         
        ++iter;
    }

    omp_set_default_device(0);
    omp_target_memcpy(u[0][0], d_u_data_0, (grid_size/2) * grid_size * grid_size * sizeof(double),
            0, 0, omp_get_initial_device(), omp_get_default_device());

    omp_set_default_device(1);
    omp_target_memcpy(u[(N/2)][0], d_u_data_1, (grid_size/2) * grid_size * grid_size * sizeof(double),
            0, 0, omp_get_initial_device(), omp_get_default_device());
    

    d_free_3d(d_u_0, d_u_data_0);
    d_free_3d(d_old_u_0, d_old_u_data_0);
    d_free_3d(d_f_0, d_f_data_0);
    d_free_3d(d_u_1, d_u_data_1);
    d_free_3d(d_old_u_1, d_old_u_data_1);
    d_free_3d(d_f_1, d_f_data_1);

    return iter;
}
