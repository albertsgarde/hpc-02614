#include<stdio.h>
#include<cblas.h>

void clearC(int m, int n, double **C) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      C[i][j] = 0;
    }
  }
}

void matmult_lib(int m, int n, int k, const double** A, const double** B, double** C) {
  clearC(m, n, C);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, *A, k, *B, n, 0.0, *C, n);
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

#define DIM1 m
#define DIM2 k
#define DIM3 n
#define IND1 i
#define IND2 l
#define IND3 j

void blkmult(int ind1_start, int ind1_end, int ind2_start, int ind2_end, int ind3_start, int ind3_end, const double** A, const double** B, double** C) {
  //printf("(%d, %d), (%d, %d), (%d, %d)\n", ind1_start, ind1_end, ind2_start, ind2_end, ind3_start, ind3_end);
  for (int IND1=ind1_start; IND1<ind1_end; ++IND1) {
    for (int IND2=ind2_start; IND2<ind2_end; ++IND2) {
      for (int IND3=ind3_start; IND3<ind3_end; ++IND3) {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
}

void matmult_blk(int m, int n, int k, const double** A, const double** B, double** C, int bs) {
  clearC(m, n, C);
  
  int num_blk1=DIM1/bs;
  int num_blk2=DIM2/bs;
  int num_blk3=DIM3/bs;

  
  for (int blk1=0; blk1<num_blk1; ++blk1) {
    for (int blk2=0; blk2<num_blk2; ++blk2) {
      for (int blk3=0; blk3<num_blk3; ++blk3) {
        blkmult(blk1*bs, blk1*bs+bs, blk2*bs, blk2*bs+bs, blk3*bs, blk3*bs+bs, A, B, C);
      }
      blkmult(blk1*bs, blk1*bs+bs, blk2*bs, blk2*bs+bs, num_blk3*bs, DIM3, A, B, C);
    }
    for (int blk3=0; blk3<num_blk3; ++blk3) {
      blkmult(blk1*bs, blk1*bs+bs, num_blk2*bs, DIM2, blk3*bs, blk3*bs+bs, A, B, C);
    }
    blkmult(blk1*bs, blk1*bs+bs, num_blk2*bs, DIM2, num_blk3*bs, DIM3, A, B, C);
  }
  for (int blk2=0; blk2<num_blk2; ++blk2) {
    for (int blk3=0; blk3<num_blk3; ++blk3) {
      blkmult(num_blk1*bs, DIM1, blk2*bs, blk2*bs+bs, blk3*bs, blk3*bs+bs, A, B, C);
    }
    blkmult(num_blk1*bs, DIM1, blk2*bs, blk2*bs+bs, num_blk3*bs, DIM3, A, B, C);
  }
  for (int blk3=0; blk3<num_blk3; ++blk3) {
    blkmult(num_blk1*bs, DIM1, num_blk2*bs, DIM2, blk3*bs, blk3*bs+bs, A, B, C);
  }
  blkmult(num_blk1*bs, DIM1, num_blk2*bs, DIM2, num_blk3*bs, DIM3, A, B, C);
}
