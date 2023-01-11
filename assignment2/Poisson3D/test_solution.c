#include <math.h>
#include "utils.h"
#include "init_u.h"

// Gives the solution for the test problem.
void test_solution(const int N, double ***u) {
    init_u(N, u, 0, true);
    for (int i = 1; i < N+1; ++i) {
        const double x = index_to_coord(N, i);
        for (int j = 1; j < N+1; ++j) {
            const double y = index_to_coord(N, j);
            for (int k = 1; k < N+1; ++k) {
                const double z = index_to_coord(N, k);
                u[i][j][k] = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
            }
        }
    }
}
