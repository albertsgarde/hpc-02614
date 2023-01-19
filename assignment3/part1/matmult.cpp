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
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, *A, k, *B, n, 0.0, *C, n);
  }
}

// STANDARD OPEN MP VERSIONS

void matmult_mkn_omp(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);

  #pragma omp parallel for
  for (int i=0; i<m; i++) {
    for (int l=0; l<k; l++) {
      for (int j=0; j<n; j++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_blk_omp(int m, int n, int k, const double** A, const double** B, double** C, int bs) {
  clearC(m, n, C);

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
}

// OFFLOADED VERSIONS

extern "C" {
  void matmult_mkn_offload(int m, int n, int k, const double** A, const double** B, double** C) {
    clearC(m, n, C);

    #pragma omp target data map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n]) 
    {
      #pragma omp target teams loop map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n]) \
                  num_teams(m) thread_limit(64)
      for (int i=0; i<m; i++) {
        #pragma omp loop bind(parallel)
        for (int l=0; l<k; l++) {
          for (int j=0; j<n; j++) {
            C[i][j] += A[i][l] * B[l][j];
          }
        }
      }
    }
  }
}

extern "C" {
  void matmult_mnk_offload(int m, int n, int k, const double** A, const double** B, double** C) {
    clearC(m, n, C);

    #pragma omp target data map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n]) 
    {
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
}

extern "C" {
  void matmult_blk_offload(int m, int n, int k, const double** A, const double** B, double** C) {
    #pragma omp target data map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n])
    {
      #pragma omp target teams loop map(to:A[:m][:k],B[:k][:n],m,n,k) map(tofrom:C[:m][:n]) \
              num_teams(2056) thread_limit(64) collapse(2)
      for (int i1=0; i1<m; i1+=BLK_GPU_BS) {
        for (int j=0; j<n; j++) {
          if (m-i1 > BLK_GPU_BS) {
            double sum[BLK_GPU_BS] = {0};
            for (int l=0; l<k; l++) {
              for (int i2=0; i2<BLK_GPU_BS; i2++) {
                sum[i2] += A[i1+i2][l] * B[l][j];
              }
              for (int i2=0; i2<BLK_GPU_BS; i2++) {
                C[i1+i2][j] = sum[i2];
              }
            } 
          } else {
            double sum[BLK_GPU_BS] = {0};
            for (int l=0; l<k; l++) {
              for (int i2=0; i2<m-i1; i2++) {
                sum[i2] += A[i1+i2][l] * B[l][j];
              }
              for (int i2=0; i2<m-i1; i2++) {
                C[i1+i2][j] = sum[i2];
              }
            }
          }
        }
      }
    }
  }
}
