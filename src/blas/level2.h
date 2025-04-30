#ifndef BLAS_LEVEL2_H
#define BLAS_LEVEL2_H
/*
This function does y := alpha*A*x + beta*y,
*/

/* 
implement blas level 2 operations here

we need:

sgemv
sgemv312
sgemv321
sgemv422
*/

__device__ __host__
void sgemv(
  bool trans, int m, int n, float alpha, const float* amat, int lda,
  const float* xvec, int incx, float beta, float* yvec, int incy
)
{
  if (!trans) {
      // y := alpha * A * x + beta * y
      for (int i = 0; i < m; ++i) {
          float temp = 0.0f;
          for (int j = 0; j < n; ++j) {
              temp += amat[i * lda + j] * xvec[j * incx];
          }
          yvec[i * incy] = alpha * temp + beta * yvec[i * incy];
      }
  } else {
      // y := alpha * A^T * x + beta * y
      for (int i = 0; i < n; ++i) {
          float temp = 0.0f;
          for (int j = 0; j < m; ++j) {
              temp += amat[j * lda + i] * xvec[j * incx];
          }
          yvec[i * incy] = alpha * temp + beta * yvec[i * incy];
      }
  }
}


void test_sgemv() {
  // Test case 1: Non-transposed matrix multiplication
  {
    const int m = 2, n = 3;
    float amat[m * n] = {1, 2, 3, 4, 5, 6};
    float xvec[n] = {1, 1, 1};
    float yvec[m] = {0, 0};
    float alpha = 1.0f, beta = 0.0f;
    int lda = n, incx = 1, incy = 1;

    sgemv(false, m, n, alpha, amat, lda, xvec, incx, beta, yvec, incy);

    assert(yvec[0] == 6.0f); // 1*1 + 2*1 + 3*1
    assert(yvec[1] == 15.0f); // 4*1 + 5*1 + 6*1
  }

  // Test case 2: Transposed matrix multiplication
  {
    const int m = 2, n = 3;
    float amat[m * n] = {1, 2, 3, 4, 5, 6};
    float xvec[m] = {1, 1};
    float yvec[n] = {0, 0, 0};
    float alpha = 1.0f, beta = 0.0f;
    int lda = n, incx = 1, incy = 1;

    sgemv(true, m, n, alpha, amat, lda, xvec, incx, beta, yvec, incy);

    assert(yvec[0] == 5.0f); // 1*1 + 4*1
    assert(yvec[1] == 7.0f); // 2*1 + 5*1
    assert(yvec[2] == 9.0f); // 3*1 + 6*1
  }

  // Test case 3: Non-transposed with beta scaling
  {
    const int m = 2, n = 3;
    float amat[m * n] = {1, 2, 3, 4, 5, 6};
    float xvec[n] = {1, 1, 1};
    float yvec[m] = {1, 1};
    float alpha = 1.0f, beta = 2.0f;
    int lda = n, incx = 1, incy = 1;

    sgemv(false, m, n, alpha, amat, lda, xvec, incx, beta, yvec, incy);

    assert(yvec[0] == 8.0f); // (1*1 + 2*1 + 3*1) + 2*1
    assert(yvec[1] == 17.0f); // (4*1 + 5*1 + 6*1) + 2*1
  }

  std::cout << "All tests passed for sgemv!" << std::endl;
}


void test_blas_level2()
{
  test_sgemv();
}

#endif