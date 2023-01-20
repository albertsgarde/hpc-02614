/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _FROBENIUS_NORM_H
#define _FROBENIUS_NORM_H

double frobenius_norm(double ***A, double ***B, const int N);

double frobenius_norm_par(double ***A, double ***B, const int N);

double frobenius_norm_gpu(double ***d_A, double ***d_B, const int N);

#endif /* _FROBENIUS_NORM_H */
