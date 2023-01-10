#include <math.h>

double frobeniusnorm(double*** U, double *** oldU, int N ){

    double norm = 0;

    for (int i = 1; i < (N + 2); i++){
        for (int j=1; j < (N + 2); j++){
            for (int k=1; k < (N + 2); k++){

                    double difference = U[i][j][k] - oldU[i][j][k];
                    norm += difference * difference
            }
        }
    }

    return Math.Sqrt(norm)

}

