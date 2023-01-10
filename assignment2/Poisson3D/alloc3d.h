#ifndef __ALLOC_3D
#define __ALLOC_3D

double ***malloc_3d(int m, int n, int k);

#define HAS_FREE_3D
void free_3d(double ***array3D, int m, int n, int k);

#endif /* __ALLOC_3D */
