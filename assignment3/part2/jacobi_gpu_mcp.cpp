/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>
#include "frobenius_norm.h"
#include "utils.h"
#include "d_alloc3d.h"

void jacobi_inner_gpu_mcp(double ***d_u, double ***d_old_u, double ***d_f, const int N) {

    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);

    #pragma omp target teams loop is_device_ptr(d_u, d_old_u, d_f)\
        collapse(3) num_teams(N*N/4) thread_limit(256)
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                double value1 = grid_spacing_sq * d_f[i][j][k];
                double value = one_sixth * (d_old_u[i-1][j][k] + d_old_u[i+1][j][k] + d_old_u[i][j-1][k] + d_old_u[i][j+1][k] 
                        + d_old_u[i][j][k-1] + d_old_u[i][j][k+1] + value1);
                d_u[i][j][k] = value;
            }
        }
    }
}

double jacobi_inner_gpu_mcp_frob(double ***d_u, double ***d_old_u, double ***d_f, const int N) {

    const double one_sixth = 1./6.;
    const double grid_spacing_sq = grid_spacing(N)*grid_spacing(N);

    double total_delta = 0;

    #pragma omp target teams loop is_device_ptr(d_u, d_old_u, d_f) reduction(+: total_delta)\
        collapse(3) num_teams(N*N/4) thread_limit(256)
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                double value1 = grid_spacing_sq * d_f[i][j][k];
                double value = one_sixth * (d_old_u[i-1][j][k] + d_old_u[i+1][j][k] + d_old_u[i][j-1][k] + d_old_u[i][j+1][k] 
                        + d_old_u[i][j][k-1] + d_old_u[i][j][k+1] + value1);
                d_u[i][j][k] = value;

                const double delta = d_u[i][j][k] - d_old_u[i][j][k];
                total_delta += delta * delta;
            }
        }
    }
    return sqrt(total_delta);
}

int jacobi_gpu_mcp(double *** u, double *** old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius) {
    int iter = 0;
    double delta_norm = INFINITY;

    double* d_u_data;
    double*** d_u = d_malloc_3d(N+2, N+2, N+2, &d_u_data);
    copy_grid_to_device(u[0][0], d_u_data, N+2);

    double* d_f_data;
    double*** d_f = d_malloc_3d(N+2, N+2, N+2, &d_f_data);
    copy_grid_to_device(f[0][0], d_f_data, N+2);
    
    double* d_old_u_data;
    double*** d_old_u = d_malloc_3d(N+2, N+2, N+2, &d_old_u_data);


    while (delta_norm > threshold && iter < iter_max) {
        double ***tmp = d_u;
        d_u = d_old_u;
        d_old_u = tmp;
        if (frobenius) {
            delta_norm = jacobi_inner_gpu_mcp_frob(d_u, d_old_u, d_f, N);
        } else {
            jacobi_inner_gpu_mcp(d_u, d_old_u, d_f, N);
        }
        ++iter;
    }

    copy_grid_from_device(u[0][0], d_u_data, N+2);

    d_free_3d(d_u, d_u_data);
    d_free_3d(d_old_u, d_old_u_data);
    d_free_3d(d_f, d_f_data);

    return iter;
}
