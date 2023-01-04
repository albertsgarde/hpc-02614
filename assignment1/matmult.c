#include<cblas.h>

void clearC(int m, int n, double **C) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      C[i][j] = 0;
    }
  }
}

void matmult_nat(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      for (int l=0; l<k; l++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_lib(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, *A, k, *B, n, 0.0, *C, n);
}

void matmult_mnk(int m, int n, int k, const double** A, const double** B, double** C) {
  matmult_nat(m, n, k, A, B, C);
}

void matmult_mkn(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  for (int i=0; i<m; i++) {
    for (int l=0; l<k; l++) {
      for (int j=0; j<n; j++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_nmk(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  for (int j=0; j<n; j++) {
    for (int i=0; i<m; i++) {
      for (int l=0; l<k; l++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_nkm(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  for (int j=0; j<n; j++) {
    for (int l=0; l<k; l++) {
      for (int i=0; i<m; i++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_kmn(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  for (int l=0; l<k; l++) {
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_knm(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  for (int l=0; l<k; l++) {
    for (int j=0; j<n; j++) {
      for (int i=0; i<m; i++) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}
