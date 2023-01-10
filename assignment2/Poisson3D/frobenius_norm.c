#include <math.h>

double frobenius_norm(double ***A, double ***B, const int N){

    double norm = 0;

    for (int i = 1; i < (N + 2); i++){
        for (int j=1; j < (N + 2); j++){
            for (int k=1; k < (N + 2); k++){

                    const double difference = A[i][j][k] - B[i][j][k];
                    norm += difference * difference;
            }
        }
    }

    return sqrt(norm);
}
