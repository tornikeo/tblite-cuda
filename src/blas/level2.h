#ifndef BLAS_LEVEL2_H
#define BLAS_LEVEL2_H
/* 
This function does y := alpha*A*x + beta*y,
*/
// template <int N>
__device__ void symv(
  int N,
  const float *amat, // 1D array representing symmetric matrix
  const float (&xvec)[],//[N],
  float (&yvec)[], //[N],
  const float beta = 0.0f,
  const float alpha = 1.0f
);

// BLAS-style matrix-vector multiplication
__device__ 
void gemv(const float* A, const float* x, float* y, 
  size_t rows, size_t cols, 
  float alpha = 1.0f, 
  float beta = 0.0f, 
  bool transpose = false);

// BLAS-style matrix-vector multiplication with a 3D matrix (first two dimensions flattened)
__device__ void gemv312(const float* A, const float* x, float* y, 
  size_t dim1, size_t dim2, size_t dim3, 
  float alpha = 1.0f, float beta = 0.0f, 
  bool transpose = false);


__device__ void gemv422(const float* A, const float* x, float* y, 
    size_t dim1, size_t dim2, size_t dim3, size_t dim4, 
    float alpha = 1.0f, float beta = 0.0f, 
    bool transpose = false);

void test_blas();
#endif