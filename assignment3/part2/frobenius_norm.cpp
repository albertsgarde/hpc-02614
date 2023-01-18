#include <math.h>

double frobenius_norm(double ***A, double ***B, const int N) {
    double norm_sq = 0;

    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                const double difference = A[i][j][k] - B[i][j][k];
                norm_sq += difference * difference;
            }
        }
    }

    return sqrt(norm_sq);
}

double frobenius_norm_par(double ***A, double ***B, const int N) {
    double norm_sq = 0;

    #pragma omp parallel for default(shared)
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                const double difference = A[i][j][k] - B[i][j][k];
                norm_sq += difference * difference;
            }
        }
    }

    return sqrt(norm_sq);
}

double frobenius_norm_gpu(double ***d_A, double ***d_B, const int N) {
    double norm_sq = 0;

    #pragma omp target teams loop is_device_ptr(d_A, d_B)\
            collapse(2) num_teams(108) thread_limit(128) reduction(+: norm_sq)
    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){
                const double difference = d_A[i][j][k] - d_B[i][j][k];
                norm_sq += difference * difference;
            }
        }
    }

    return sqrt(norm_sq);
}
