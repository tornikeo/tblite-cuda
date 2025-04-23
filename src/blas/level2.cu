#include "../limits.h"

/* TODO: priority low Consider this. Should this be a templated function */
__device__ void symv(
  int N,
  const float *amat, // 1D array representing symmetric matrix
  const float (&xvec)[],//[N],
  float (&yvec)[], //[N],
  const float beta = 0.0f,
  const float alpha = 1.0f
)
{
  for (int i = 0; i < N; i++) {
   float temp = 0.0f;
   for (int j = 0; j < N; j++) {
    temp += amat[i * N + j] * xvec[j]; // Map 1D array to 2D matrix
   }
   yvec[i] = alpha * temp + beta * yvec[i];
  }
}
__device__ 
void gemv(const float* A, const float* x, float* y, 
  size_t rows, size_t cols, 
  float alpha = 1.0f, 
  float beta = 0.0f, 
  bool transpose = false) {
    size_t row, col;
    float result;

    if (!transpose) {
        // y := alpha * A * x + beta * y
        for (row = 0; row < rows; ++row) {
            result = 0.0f;
            for (col = 0; col < cols; ++col) {
                result += A[row * cols + col] * x[col];
            }
            y[row] = alpha * result + beta * y[row];
        }
    } else {
        // y := alpha * A^T * x + beta * y
        for (col = 0; col < cols; ++col) {
            result = 0.0f;
            for (row = 0; row < rows; ++row) {
                result += A[row * cols + col] * x[row];
            }
            y[col] = alpha * result + beta * y[col];
        }
    }
}

// BLAS-style matrix-vector multiplication with a 3D matrix (first two dimensions flattened)
__device__ void gemv312(const float* A, const float* x, float* y, 
  size_t dim1, size_t dim2, size_t dim3, 
  float alpha = 1.0f, float beta = 0.0f, 
  bool transpose = false) {
    size_t flattened_rows = dim1 * dim2;
    size_t cols = dim3;

    if (!transpose) {
        // Flatten the first two dimensions and call gemv
        gemv(A, x, y, flattened_rows, cols, alpha, beta, false);
    } else {
        // Transpose case: treat A as if it were transposed
        gemv(A, x, y, cols, flattened_rows, alpha, beta, true);
    }
}


/* */