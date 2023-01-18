/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */


int jacobi_gpu_mcp(double ***u, double ***old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius);
int jacobi_gpu_mcp_split_frob(double ***u, double ***old_u, double ***f, const int N, const int iter_max, const double threshold, const bool frobenius);
