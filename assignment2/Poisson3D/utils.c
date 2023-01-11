#include <math.h>

int coord_to_index(const int N, const double coord) {
    return (int) (((int)N+1.)*(coord + 1.)/2.);
}

double index_to_coord(const int N, const int index) {
    return 2.*((double)index)/((double)N+1.) - 1.;
}

double grid_spacing(const int N) {
    return 2./(double)(N+1);
}

void init_edges(const int N, double ***a) {
    for (int i = 0; i < N+2; ++i) {
        a[i][0][0] = NAN;
        a[i][0][N+1] = NAN;
        a[i][N+1][0] = NAN;
        a[i][N+1][N+1] = NAN;
    }
    for (int j = 0; j < N+2; ++j) {
        a[0][j][0] = NAN;
        a[0][j][N+1] = NAN;
        a[N+1][j][0] = NAN;
        a[N+1][j][N+1] = NAN;
    }
    for (int k = 0; k < N+2; ++k) {
        a[0][0][k] = NAN;
        a[0][N+1][k] = NAN;
        a[N+1][0][k] = NAN;
        a[N+1][N+1][k] = NAN;
    }
}

void subtract_arrays(const int N, double ***A, double ***B, double ***C) {
    for (int i = 0; i < N+2; ++i) {
        for (int j = 0; j < N+2; ++j) {
            for (int k = 0; k < N+2; ++k) {
                C[i][j][k] = A[i][j][k] - B[i][j][k];
            }
        }
    }
}
