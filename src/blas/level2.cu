#include "../limits.h"
#include "../utils/gpu.h"
#include <stdio.h>
#include <iostream>
#include <cassert>
#include "level2.h"

// #ifndef TEST_LEVEL2_H
// #define TEST_LEVEL2_H

void test_blas() {
  // Test case for gemv with 2D matrix
  {
    constexpr int N = 2, M = 3;
    float A[N][M] = {{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}};
    float x[M] = {1.0f, 2.0f, 3.0f};
    float y[N] = {0.0f, 0.0f};
    float expected_y[N] = {14.0f, 32.0f}; // Computed manually: A*x

    gemv<N, M>(A, x, y);

    for (int i = 0; i < N; ++i) {
      assert(fabs(y[i] - expected_y[i]) < 1e-5);
    }
    std::cout << "Test case 1 (2D gemv) passed!" << std::endl;
  }

  // Test case for gemv with 3D matrix
  {
    constexpr int N = 2, K = 2, M = 2;
    float A[N][K][M] = {{{1.0f, 2.0f}, {3.0f, 4.0f}}, {{5.0f, 6.0f}, {7.0f, 8.0f}}};
    float x[K][M] = {{1.0f, 2.0f}, {3.0f, 4.0f}};
    float y[N] = {0.0f, 0.0f};
    float expected_y[N] = {50.0f, 114.0f}; // Computed manually: A*x

    gemv<N, K, M>(A, x, y);

    for (int i = 0; i < N; ++i) {
      assert(fabs(y[i] - expected_y[i]) < 1e-5);
    }
    std::cout << "Test case 2 (3D gemv) passed!" << std::endl;
  }

  // Test case for gemv with 4D matrix
  {
    constexpr int M = 2, N = 2;
    float A[M][N][M][N] = {
      {{{1.0f, 2.0f}, {3.0f, 4.0f}}, {{5.0f, 6.0f}, {7.0f, 8.0f}}},
      {{{9.0f, 10.0f}, {11.0f, 12.0f}}, {{13.0f, 14.0f}, {15.0f, 16.0f}}}
    };
    float x[M][N] = {{1.0f, 2.0f}, {3.0f, 4.0f}};
    float y[M][N] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
    float expected_y[M][N] = {{250.0f, 260.0f}, {618.0f, 644.0f}}; // Computed manually: A*x

    gemv<M, N>(A, x, y);

    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        assert(fabs(y[i][j] - expected_y[i][j]) < 1e-5);
      }
    }
    std::cout << "Test case 3 (4D gemv) passed!" << std::endl;
  }
}


// /* TODO: priority low Consider this. Should this be a templated function */
// __device__ void symv(
//     int N,
//     const float *amat,     // 1D array representing symmetric matrix
//     const float (&xvec)[], //[N],
//     float (&yvec)[],       //[N],
//     const float beta = 0.0f,
//     const float alpha = 1.0f)
// {
//   for (int i = 0; i < N; i++)
//   {
//     float temp = 0.0f;
//     for (int j = 0; j < N; j++)
//     {
//       temp += amat[i * N + j] * xvec[j]; // Map 1D array to 2D matrix
//     }
//     yvec[i] = alpha * temp + beta * yvec[i];
//   }
// }

// __device__ void gemv(const float *A, const float *x, float *y,
//                      size_t rows, size_t cols,
//                      float alpha = 1.0f,
//                      float beta = 0.0f,
//                      bool transpose = false)
// {
//   size_t row, col;
//   float result;

//   if (!transpose)
//   {
//     // y := alpha * A * x + beta * y
//     for (row = 0; row < rows; ++row)
//     {
//       result = 0.0f;
//       for (col = 0; col < cols; ++col)
//       {
//         result += A[row * cols + col] * x[col];
//       }
//       y[row] = alpha * result + beta * y[row];
//     }
//   }
//   else
//   {
//     // y := alpha * A^T * x + beta * y
//     for (col = 0; col < cols; ++col)
//     {
//       result = 0.0f;
//       for (row = 0; row < rows; ++row)
//       {
//         result += A[row * cols + col] * x[row];
//       }
//       y[col] = alpha * result + beta * y[col];
//     }
//   }
// }

// __global__ void test_gemv()
// {
//   // Define a small matrix A (2x3), vector x (3), and vector y (2)
//   const int rows = 2, cols = 3;
//   float A[rows * cols] = {1.0f, 2.0f, 3.0f,
//                           4.0f, 5.0f, 6.0f};
//   float x[cols] = {1.0f, 1.0f, 1.0f};
//   float y[rows] = {0.0f, 0.0f};

//   // Call the __device__ gemv function directly
//   gemv(A, x, y, rows, cols);

//   // Print the result
//   printf("Result y: [%f, %f]\n", y[0], y[1]);

//   // Expected results
//   float expected_y[rows] = {6.0f, 15.0f};

//   // Assert the results
//   assert(fabs(y[0] - expected_y[0]) < 1e-5 && "Assertion failed for y[0]");
//   assert(fabs(y[1] - expected_y[1]) < 1e-5 && "Assertion failed for y[1]");
// }

// // BLAS-style matrix-vector multiplication with a 3D matrix (first two dimensions flattened)
// __device__ void gemv312(const float *A, const float *x, float *y,
//                         size_t dim1, size_t dim2, size_t dim3,
//                         float alpha = 1.0f, float beta = 0.0f,
//                         bool transpose = false)
// {
//   size_t flattened_rows = dim1 * dim2;
//   size_t cols = dim3;
//   // printf("gemv312 called with A(%i %i %i) @ x(%i) + y(%i %i)\n", dim1, dim2, dim3, )
//   if (!transpose)
//   {
//     // Flatten the first two dimensions and call gemv
//     gemv(A, x, y, flattened_rows, cols, alpha, beta, false);
//   }
//   else
//   {
//     // Transpose case: treat A as if it were transposed
//     gemv(A, x, y, cols, flattened_rows, alpha, beta, true);
//   }
// }

// __global__ void test_gemv312()
// {
//   // Define a small 3D matrix A (2x2x3), vector x (3), and vector y (4)
//   const int dim1 = 2, dim2 = 2, dim3 = 3;
//   const int flattened_rows = dim1 * dim2;
//   float A[flattened_rows * dim3] = {
//       1.0f, 2.0f, 3.0f,
//       4.0f, 5.0f, 6.0f,
//       7.0f, 8.0f, 9.0f,
//       10.0f, 11.0f, 12.0f};
//   float x[dim3] = {1.0f, 1.0f, 1.0f};
//   float y[flattened_rows] = {0.0f, 0.0f, 0.0f, 0.0f};

//   // Call the __device__ gemv312 function directly
//   gemv312(A, x, y, dim1, dim2, dim3);

//   // Print the result
//   printf("Result y: [%f, %f, %f, %f]\n", y[0], y[1], y[2], y[3]);

//   // Expected results
//   float expected_y[flattened_rows] = {6.0f, 15.0f, 24.0f, 33.0f};

//   // Assert the results
//   for (int i = 0; i < flattened_rows; ++i)
//   {
//     assert(fabs(y[i] - expected_y[i]) < 1e-5 && "Assertion failed for y[i]");
//   }
// }

// /*
// subroutine wrap_dgemv321(amat, xvec, yvec, alpha, beta, trans)
//    real(dp), intent(in), contiguous, target :: amat(:, :, :)
//    real(dp), intent(in), contiguous, target :: xvec(:, :)
//    real(dp), intent(inout) :: yvec(:)
//    real(dp), intent(in), optional :: alpha
//    real(dp), intent(in), optional :: beta
//    character(len=1), intent(in), optional :: trans
//    real(dp), pointer :: aptr(:, :), xptr(:)
//    character(len=1) :: tra
//    if (present(trans)) then
//       tra = trans
//    else
//       tra = 'n'
//    end if
//    if (any(tra == ['n', 'N'])) then
//       aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
//       xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
//    else
//       aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
//       xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
//    end if
//    call wrap_gemv(aptr, xptr, yvec, alpha, beta, tra)
// end subroutine wrap_dgemv321
// */
// __device__ void gemv321(const float *A, const float *x, float *y,
//   size_t dim1, size_t dim2, size_t dim3,
//   float alpha = 1.0f, float beta = 0.0f,
//   bool transpose = false)
// {
//   size_t rows = dim1;
//   size_t flattened_cols = dim2 * dim3;

//   if (!transpose)
//   {
//     // Flatten the last two dimensions and call gemv
//     gemv(A, x, y, rows, flattened_cols, alpha, beta, false);
//   }
//   else
//   {
//     // Transpose case: treat A as if it were transposed
//     gemv(A, x, y, flattened_cols, rows, alpha, beta, true);
//   }
// }

// /*
// subroutine wrap_dgemv422(amat, xvec, yvec, alpha, beta, trans)
//    real(dp), intent(in), contiguous, target :: amat(:, :, :, :)
//    real(dp), intent(in), contiguous, target :: xvec(:, :)
//    real(dp), intent(inout), contiguous, target :: yvec(:, :)
//    real(dp), intent(in), optional :: alpha
//    real(dp), intent(in), optional :: beta
//    character(len=1), intent(in), optional :: trans
//    real(dp), pointer :: aptr(:, :), yptr(:), xptr(:)
//    character(len=1) :: tra
//    if (present(trans)) then
//       tra = trans
//    else
//       tra = 'n'
//    end if
//    aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)*size(amat, 4)) => amat
//    xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
//    yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
//    call wrap_gemv(aptr, xptr, yptr, alpha, beta, tra)
// end subroutine wrap_dgemv422
// */

// __device__ void gemv422(const float *A, const float *x, float *y,
//                         size_t dim1, size_t dim2, size_t dim3, size_t dim4,
//                         float alpha = 1.0f, float beta = 0.0f,
//                         bool transpose = false)
// {
//   size_t flattened_rows = dim1 * dim2;
//   size_t flattened_cols = dim3 * dim4;

//   if (!transpose)
//   {
//     // Flatten the first two dimensions and the last two dimensions, then call gemv
//     gemv(A, x, y, flattened_rows, flattened_cols, alpha, beta, false);
//   }
//   else
//   {
//     // Transpose case: treat A as if it were transposed
//     gemv(A, x, y, flattened_cols, flattened_rows, alpha, beta, true);
//   }
// }

// __global__ void test_gemv422()
// {
//   // Define a small 4D matrix A (2x2x2x3), vector x (6), and vector y (4)
//   const int dim1 = 2, dim2 = 2, dim3 = 2, dim4 = 3;
//   const int flattened_rows = dim1 * dim2;
//   const int flattened_cols = dim3 * dim4;
//   float A[flattened_rows * flattened_cols] = {
//       1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f,
//       7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f,
//       13.0f, 14.0f, 15.0f, 16.0f, 17.0f, 18.0f,
//       19.0f, 20.0f, 21.0f, 22.0f, 23.0f, 24.0f};
//   float x[flattened_cols] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
//   float y[flattened_rows] = {0.0f, 0.0f, 0.0f, 0.0f};

//   // Call the __device__ gemv422 function directly
//   gemv422(A, x, y, dim1, dim2, dim3, dim4);

//   // Print the result
//   printf("Result y: [%f, %f, %f, %f]\n", y[0], y[1], y[2], y[3]);

//   // Expected results
//   float expected_y[flattened_rows] = {21.0f, 57.0f, 93.0f, 129.0f};

//   // Assert the results
//   for (int i = 0; i < flattened_rows; ++i)
//   {
//     assert(fabs(y[i] - expected_y[i]) < 1e-5 && "Assertion failed for y[i]");
//   }
// }

// void test_blas()
// {
//   cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
//   cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
//   test_gemv<<<1, 1>>>();
//   test_gemv312<<<1, 1>>>();
//   test_gemv422<<<1, 1>>>();
//   gpuErrchk(cudaPeekAtLastError());
//   gpuErrchk(cudaDeviceSynchronize());
// }
