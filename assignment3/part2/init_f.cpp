#include <math.h>
#include <stdbool.h>
#include "utils.h"

void init_f_ass(const int N, double ***f) {
    for (int i = 1; i < N+1; ++i) {
        for (int j = 1; j < N+1; ++j) {
            for (int k = 1; k < N+1; ++k) {
                f[i][j][k] = 0.;
            }
        }
    }

    const double max_radiator_x = -3./8.;
    const double max_radiator_y = -1./2.;
    const double min_radiator_z = -2./3.;
    const double max_radiator_z = 0.;

    const int max_radiator_i = coord_to_index(N, max_radiator_x);
    const int max_radiator_j = coord_to_index(N, max_radiator_y);
    const int min_radiator_k = coord_to_index(N, min_radiator_z);
    const int max_radiator_k = coord_to_index(N, max_radiator_z);

    for (int i = 0; i < max_radiator_i; ++i) {
        for(int j = 0; j < max_radiator_j; ++j) {
            for (int k = min_radiator_k; k < max_radiator_k; ++k) {
                f[i][j][k] = 200.;
            }
        }
    }
}

void init_f_test(const int N, double ***f) {
    for (int i = 1; i < N+1; ++i) {
        for (int j = 1; j < N+1; ++j) {
            for (int k = 1; k < N+1; ++k) {
                double x = index_to_coord(N, i);
                double y = index_to_coord(N, j);
                double z = index_to_coord(N, k);
                f[i][j][k] = 3*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
            }
        }
    }
}

void init_f(const int N, double ***f, const bool test) {
    if (test) {
        init_f_test(N, f);
    } else {
        init_f_ass(N, f);
    }
    // Set edges to NAN, since they should never be used.
    init_edges(N, f);
}
