/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H

void
gauss_seidel(double *const *const *const U, const double *const *const *const f, const int N, const int iter_max, const double threshold);

#endif
