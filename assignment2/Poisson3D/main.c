/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc3d.h"
#include "print.h"

#ifdef _JACOBI
#include "jacobi.h"
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#include "init_u.h"
#include "init_f.h"
#include "utils.h"
#include "test_solution.h"
#include "frobenius_norm.h"

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

    char *output_prefix = "poisson_res";
    char *solution_output_prefix = "poisson_sol";
    char *error_output_prefix = "poisson_err";
    char *output_ext    = "";


    /* get the paramters from the command line */
    if (argc < 5) {
        printf("Usage: %s N iter_max tolerance start_T [output_type]\n", argv[0]);
        return(1);
    }
    const int N         = atoi(argv[1]);	// grid size
    const int iter_max  = atoi(argv[2]);  // max. no. of iterations
    const double tolerance = atof(argv[3]);  // tolerance
    const double start_T   = atof(argv[4]);  // start T for all inner grid points
    // Whether to use test boundary conditions and f.
    bool test = false;
    int output_type = 0;
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-t") == 0) {
            test = true;
        } else {
            output_type = atoi(argv[i]);
        }
    }

    // allocate memory
    double *** u = NULL;
    if ( (u = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }

    double *** f = NULL;
    if ( (f = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array f: allocation failed");
        exit(-1);
    }

    // initialize grid with boundary conditions

    init_u(N, u, start_T, test);
    init_f(N, f, test);

    
    #ifdef _JACOBI
    double *** old_u = NULL;
    if ( (old_u = malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array old_u: allocation failed");
        exit(-1);
    }
    init_u(N, old_u, start_T, test);
    jacobi(u, old_u, f, N, iter_max, tolerance);
    #endif

    #ifdef _GAUSS_SEIDEL
    gauss_seidel(u, f, N, iter_max, tolerance);
    #endif

    // dump  results if wanted 

    char output_filename[FILENAME_MAX];
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: \n", output_filename);
	    print_binary(output_filename, N+2, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: \n", output_filename);
	    print_vtk(output_filename, N+2, u);

        if (test) {
            double *** solution = NULL;
            if ( (solution = malloc_3d(N+2, N+2, N+2)) == NULL ) {
                perror("array solution: allocation failed");
                exit(-1);
            }
            test_solution(N, solution);
            char solution_output_filename[FILENAME_MAX];
	        sprintf(solution_output_filename, "%s_%d%s", solution_output_prefix, N, output_ext);
            fprintf(stderr, "Write solution VTK file to %s: \n", solution_output_filename);
            print_vtk(solution_output_filename, N+2, solution);

            double *** error_array = NULL;
            if ( (error_array = malloc_3d(N+2, N+2, N+2)) == NULL ) {
                perror("array error_array: allocation failed");
                exit(-1);
            }
            subtract_arrays(N, u, solution, error_array);
            char error_output_filename[FILENAME_MAX];
	        sprintf(error_output_filename, "%s_%d%s", error_output_prefix, N, output_ext);
            fprintf(stderr, "Write error VTK file to %s: \n", error_output_filename);
            print_vtk(error_output_filename, N+2, error_array);

            const double error = frobenius_norm(u, solution, N);
            printf("Error: %f\n", error);
        }
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free_3d(u);

    return(0);
}
