#include <stdlib.h>
#include <omp.h>

double*** d_malloc_3d(int m, int n, int k, double** data) {
    if (m <= 0 || n <= 0 || k <= 0)
        return NULL;

    double ***p = (double***)omp_target_alloc(m * sizeof(double **) + m * n * sizeof(double *), omp_get_default_device());

    if (p == NULL)
        return NULL;

    double *a = (double*)omp_target_alloc(m * n * k * sizeof(double), omp_get_default_device());
    if (a == NULL) {
        omp_target_free(p, omp_get_default_device());
        return NULL;
    }
    #pragma omp target is_device_ptr(p, a)
    for (int i = 0; i < m; i++)
        p[i] = (double **) p + m + i * n;
    
    #pragma omp target is_device_ptr(p, a)
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            p[i][j] = a + (i * n * k) + (j * k);

    *data = a;
    return p;
}

void d_free_3d(double*** p, double* data) {
    omp_target_free(data,
            omp_get_default_device());
    omp_target_free(p,
            omp_get_default_device());
}
