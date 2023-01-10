#include <math.h>
#include <stdbool.h>

void init_u_internal(const int N, double ***const u, const double start_T) {
    for (int i = 1; i < N+1; ++i) {
        for (int j = 1; j < N+1; ++j) {
            for (int k = 1; k < N+1; ++k) {
                u[i][j][k] = start_T;
            }
        }
    }
}

void init_u_corners(const int N, double ***const u) {
    for (int i = 0; i < N+2; ++i) {
        u[i][0][0] = NAN;
        u[i][0][N+1] = NAN;
        u[i][N+1][0] = NAN;
        u[i][N+1][N+1] = NAN;
    }
    for (int j = 0; j < N+2; ++j) {
        u[0][j][0] = NAN;
        u[0][j][N+1] = NAN;
        u[N+1][j][0] = NAN;
        u[N+1][j][N+1] = NAN;
    }
    for (int k = 0; k < N+2; ++k) {
        u[0][0][k] = NAN;
        u[0][N+1][k] = NAN;
        u[N+1][0][k] = NAN;
        u[N+1][N+1][k] = NAN;
    }
}

void init_u_boundary_ass(const int N, double ***const u) {
    for (int i = 1; i < N+1; ++i) {
        for (int j = 1; j < N+1; ++j) {
            u[i][j][0] = 20.;
            u[i][j][N+1] = 20.;
        }
        for (int k = 1; k < N+1; ++k) {
            u[i][0][k] = 0.;
            u[i][N+1][k] = 20.;
        }
    }
    for (int j = 1; j < N+1; ++j) {
        for (int k = 1; k < N+1; ++k) {
            u[0][j][k] = 20.;
            u[N+1][j][k] = 20.;
        }
    }
}

void init_u_boundary_test(const int N, double ***const u) {
    for (int i = 1; i < N+1; ++i) {
        for (int j = 1; j < N+1; ++j) {
            u[i][j][0] = 0.;
            u[i][j][N+1] = 0.;
        }
        for (int k = 1; k < N+1; ++k) {
            u[i][0][k] = 0.;
            u[i][N+1][k] = 0.;
        }
    }
    for (int j = 1; j < N+1; ++j) {
        for (int k = 1; k < N+1; ++k) {
            u[0][j][k] = 0.;
            u[N+1][j][k] = 0.;
        }
    }
}

int init_u(const int N, double ***const u, const double start_T, const bool test) {
    init_u_internal(N, u, start_T);
    init_u_corners(N, u);
    if (test) {
        init_u_boundary_test(N, u);
    } else {
        init_u_boundary_ass(N, u);
    }
}
