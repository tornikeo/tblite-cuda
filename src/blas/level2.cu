#include "../limits.h"
#include <stdio.h>
#include <cassert>

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

__global__ void test_gemv() {
  // Define a small matrix A (2x3), vector x (3), and vector y (2)
  const int rows = 2, cols = 3;
  float A[rows * cols] = {1.0f, 2.0f, 3.0f, 
                          4.0f, 5.0f, 6.0f};
  float x[cols] = {1.0f, 1.0f, 1.0f};
  float y[rows] = {0.0f, 0.0f};

  // Call the __device__ gemv function directly
  gemv(A, x, y, rows, cols);

  // Print the result
  printf("Result y: [%f, %f]\n", y[0], y[1]);

  // Expected results
  float expected_y[rows] = {6.0f, 15.0f};

  // Assert the results
  assert(fabs(y[0] - expected_y[0]) < 1e-5 && "Assertion failed for y[0]");
  assert(fabs(y[1] - expected_y[1]) < 1e-5 && "Assertion failed for y[1]");
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

__global__ void test_gemv312() {
    // Define a small 3D matrix A (2x2x3), vector x (3), and vector y (4)
    const int dim1 = 2, dim2 = 2, dim3 = 3;
    const int flattened_rows = dim1 * dim2;
    float A[flattened_rows * dim3] = {
        1.0f, 2.0f, 3.0f, 
        4.0f, 5.0f, 6.0f, 
        7.0f, 8.0f, 9.0f, 
        10.0f, 11.0f, 12.0f
    };
    float x[dim3] = {1.0f, 1.0f, 1.0f};
    float y[flattened_rows] = {0.0f, 0.0f, 0.0f, 0.0f};

    // Call the __device__ gemv312 function directly
    gemv312(A, x, y, dim1, dim2, dim3);

    // Print the result
    printf("Result y: [%f, %f, %f, %f]\n", y[0], y[1], y[2], y[3]);

    // Expected results
    float expected_y[flattened_rows] = {6.0f, 15.0f, 24.0f, 33.0f};

    // Assert the results
    for (int i = 0; i < flattened_rows; ++i) {
        assert(fabs(y[i] - expected_y[i]) < 1e-5 && "Assertion failed for y[i]");
    }
}

void test_blas() {
    test_gemv<<<1, 1>>>();
    test_gemv312<<<1, 1>>>();
}
