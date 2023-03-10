#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<cublas_v2.h>
#include<cuda_runtime_api.h>

extern "C" { 
  #include<cblas.h> 
}

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define BLK_GPU_BS 8

void clearC(int m, int n, double **C) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      C[i][j] = 0;
    }
  }
}

extern "C" {
  void matmult_lib(int m, int n, int k, const double** A, const double** B, double** C) {
    double t = omp_get_wtime();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, *A, k, *B, n, 0.0, *C, n);
    t = omp_get_wtime() - t;
    printf("elapsed = %3.5f\n", t);
  }
}

// STANDARD OPEN MP VERSIONS

extern "C" {
  void matmult_mkn_omp(int m, int n, int k, const double** A, const double** B, double** C) {
    clearC(m, n, C);

    double t = omp_get_wtime();
    #pragma omp parallel for
    for (int i=0; i<m; i++) {
      for (int l=0; l<k; l++) {
        for (int j=0; j<n; j++) {
          C[i][j] += A[i][l] * B[l][j];
        }
      }
    }
    t = omp_get_wtime() - t;
    printf("elapsed = %3.5f\n", t);
  }
}

extern "C" {
void matmult_blk_omp(int m, int n, int k, const double** A, const double** B, double** C, int bs) {
  clearC(m, n, C);

  double t = omp_get_wtime();
  #pragma omp parallel for collapse(3)
  for (int i1=0; i1<m; i1+=bs) {
    for (int l1=0; l1<k; l1+=bs) {
      for (int j1=0; j1<n; j1+=bs) {
        for (int i2=0; i2<MIN(m-i1, bs); i2++) {
          for (int l2=0; l2<MIN(k-l1, bs); l2++) {
            for (int j2=0; j2<MIN(n-j1, bs); j2++) {
              C[i1+i2][j1+j2] += A[i1+i2][l1+l2] * B[l1+l2][j1+j2];
            }
          }
        }
      }
    }
  }
  t = omp_get_wtime() - t;
  printf("elapsed = %3.5f\n", t);
}
}

// OFFLOADED VERSIONS

extern "C" {
void matmult_mkn_offload(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);

  #pragma omp target teams loop map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n]) \
          num_teams(m)
  for (int i=0; i<m; i++) {
    for (int l=0; l<k; l++) {
      for (int j=0; j<n; j++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}
}

extern "C" {
void matmult_mnk_offload(int m, int n, int k, const double** A, const double** B, double** C) {
  #pragma omp target teams loop map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n]) \
                num_teams(m) thread_limit(64)
  for (int i=0; i<m; i++) {
    #pragma omp loop bind(parallel)
    for (int j=0; j<n; j++) {
      double sum = 0;
      for (int l=0; l<k; l++) {
        sum += A[i][l] * B[l][j];
      }
      C[i][j] = sum;
    }
  }
}
}

extern "C" {
  void matmult_blk_offload(int m, int n, int k, const double** A, const double** B, double** C) { 
    double t = omp_get_wtime();
    #pragma omp target data map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n])
    {  
    double t2 = omp_get_wtime();
      #pragma omp target teams loop map(to:A[:m][:k],B[:k][:n],m,n,k) \
              map(tofrom:C[:m][:n]) collapse(2)
      for (int i1=0; i1<m; i1+=BLK_GPU_BS) {
        for (int j=0; j<n; j++) {
          if (BLK_GPU_BS < m-i1) {
            double sum[BLK_GPU_BS] = {0};
            for (int l=0; l<k; l++) {
              for (int i2=0; i2<BLK_GPU_BS; i2++) {
                sum[i2] += A[i1+i2][l] * B[l][j];
              }
            }
            for (int i2=0; i2<BLK_GPU_BS; i2++) {
              C[i1+i2][j] = sum[i2];
            }
          } else {
            for (int i2=0; i2<m-i1; i2++) {
              double sum = 0.0;
              for (int l=0; l<k; l++) {
                sum += A[i1+i2][l] * B[l][j];
              }
              C[i1+i2][j] = sum;
            }
          }
        }
      }
      printf("time w/o transfer = %3.5f seconds\n", omp_get_wtime() - t2);
    }
    printf("time w/ transfer = %3.5f seconds\n", omp_get_wtime() - t);
  }
}

// ASYNC OFFLOAD VERSION

extern "C" {
  void matmult_blk_offload_inner(int m, int n, int k, const double** A, const double** B, double** C, const int num_teams) {    
    #pragma omp target teams loop nowait map(to:A[:m][:k],B[:k][:n],m,n,k) \
            map(tofrom:C[:m][:n]) collapse(2) num_teams(m/BLK_GPU_BS*n) thread_limit(64)
    for (int i1=0; i1<m; i1+=BLK_GPU_BS) {
      for (int j=0; j<n; j++) {
        if (BLK_GPU_BS < m-i1) {
          double sum[BLK_GPU_BS] = {0};
          for (int l=0; l<k; l++) {
            for (int i2=0; i2<BLK_GPU_BS; i2++) {
              sum[i2] += A[i1+i2][l] * B[l][j];
            }
          }
          for (int i2=0; i2<BLK_GPU_BS; i2++) {
            C[i1+i2][j] = sum[i2];
          }
        } else {
          for (int i2=0; i2<m-i1; i2++) {
            double sum = 0.0;
            for (int l=0; l<k; l++) {
              sum += A[i1+i2][l] * B[l][j];
            }
            C[i1+i2][j] = sum;
          }
        }
      }
    }
  }
}

#define NUM_SLABS 8

extern "C" {
  void matmult_asy_offload(int m, int n, int k, const double** A, const double** B, double** C) {
    #pragma omp target data map(to:B[:k][:n])
    {
    for (int slab = 0; slab < NUM_SLABS; slab++) {
      const int slab_start = slab * m / NUM_SLABS;
      const int slab_length = (slab+1)* m / NUM_SLABS-slab_start;

      const double** a = &A[slab_start];
      double** c = &C[slab_start];

      matmult_blk_offload_inner(slab_length, n, k, a, B, c, m);
    }
    }
    #pragma omp taskwait
  }
}

// LIB OFFLOAD VERSION

extern "C" {
  void matmult_lib_offload(int m, int n, int k, const double** A, const double** B, double** C) {
    cublasHandle_t handle;
    cublasCreate(&handle);
    double alpha = 1.0;
    double beta = 0.0;
    const double* a = *A;
    const double* b = *B;
    double* c = *C;
    #pragma omp target data map(to:a[:m*k],b[:k*n]) map(tofrom:c[:m*n]) use_device_ptr(a,b,c)
    {
      cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, a, k, b, n, &beta, c, n);
    }
    cublasDestroy(handle);
  }
}