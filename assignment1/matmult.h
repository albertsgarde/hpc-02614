void matmult_nat(int m, int n, int k, const double** A, const double** B, double** C);
void matmult_lib(int m, int n, int k, const double** A, const double** B, double** C);

void matmult_mnk(int m, int n, int k, const double** A, const double** B, double** C);
void matmult_mkn(int m, int n, int k, const double** A, const double** B, double** C);
void matmult_nmk(int m, int n, int k, const double** A, const double** B, double** C);
void matmult_nkm(int m, int n, int k, const double** A, const double** B, double** C);
void matmult_kmn(int m, int n, int k, const double** A, const double** B, double** C);
void matmult_knm(int m, int n, int k, const double** A, const double** B, double** C);
