/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBI_H
#define _JACOBI_H

int jacobi(double *const *const *const u, double *const *const *const old_u, const double *const *const *const f, const int N, const int iter_max, const double threshold);

#endif
