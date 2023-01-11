#include <math.h>

double frobenius_norm(double ***A, double ***B, const int N){
    double norm = 0;

    for (int i = 1; i < (N + 1); i++){
        for (int j=1; j < (N + 1); j++){
            for (int k=1; k < (N + 1); k++){

                    const double difference = A[i][j][k] - B[i][j][k];
                    norm += difference * difference;
            }
        }
    }

    return sqrt(norm);
}
