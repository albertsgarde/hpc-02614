/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */


int jacobi_gpu_async(double ***u, double ***old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius);
int jacobi_gpu_async_alt(double ***u, double ***old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius);
