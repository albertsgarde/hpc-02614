/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H

int
gauss_seidel(double ***u, double ***f, const int N, const int iter_max, const double threshold);

#endif
