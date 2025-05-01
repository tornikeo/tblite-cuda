#include "../limits.h"
#include "../utils/gpu.h"
#include "../utils/array.h"
#include <stdio.h>
#include <cassert>
#include <curand.h>
#include <cmath>


__device__ __host__
void symv(
  const float (&amat)[MAX_NSH][MAX_NSH],
  const float (&xvec)[MAX_NSH],
  float (&yvec)[MAX_NSH],
  const float alpha,
  const float beta
)
{
  for (int i = 0; i < MAX_NSH; i++)
  {
    float temp = 0.0f;
    for (int j = 0; j < MAX_NSH; j++)
    {
      temp += amat[i][j] * xvec[j];
    }
    yvec[i] = alpha * temp + beta * yvec[i];
  }
}

template <int N, int M, int L, int K>
__device__ __host__
void gemv( 
  const float (&amat)[N][M],
  const float (&xvec)[L],
  float (&yvec)[K],
  const float alpha,
  const float beta,
  const bool trans
)
{
  if(trans) // M,N @ L + K
  {
    assert(N == L);
    assert(M == K);
    for (int i = 0; i < M; i++) {
      float temp = 0.0f;
      for (int j = 0; j < N; j++) {
        temp += amat[j][i] * xvec[j];
      }
      yvec[i] = alpha * temp + beta * yvec[i];
    }
  } else { // N,M @ L + K
    assert(M == L);
    assert(N == K);
    for (int i = 0; i < N; i++) {
      float temp = 0.0f;
      for (int j = 0; j < M; j++) {
        temp += amat[i][j] * xvec[j];
      }
      yvec[i] = alpha * temp + beta * yvec[i];
    }
  }
}

__global__ 
void test_gemv_0()
{
  printf("========================================\n");
  float amat[3][2] {0}; arange(amat);
  float xvec[3] {0, 1, 1};
  float yvec[2] {0};

  const float alpha = 1.0f;
  const float beta = 1.0f;
  const bool trans = true;
  printf("\n amat = "); printr(amat);
  printf("\n xvec = "); printr(xvec); 
  gemv(amat, xvec, yvec, alpha, beta, trans);
  printf("========== result =========");

  // Expected results
  float expected_yvec[2] {6, 8};
  printf("\n result = "); printr(yvec);

  for (int i = 0; i < 2; ++i)
  {
    assert(fabs(yvec[i] - expected_yvec[i]) < 1e-5 && "Assertion failed for yvec[i]");
  }
}

__global__
void test_gemv_1()
{
  printf("========================================\n");
  float amat[3][2] {0}; arange(amat);
  float xvec[2] {1, 1};
  float yvec[3] {0};

  const float alpha = 1.0f;
  const float beta = 1.0f;
  const bool trans = false;

  printf("\n amat = "); printr(amat);
  printf("\n xvec = "); printr(xvec); 
  gemv(amat, xvec, yvec, alpha, beta, trans);
  printf("========== result =========");

  // Expected results
  float expected_yvec[3] {1, 5, 9};
  printf("\n result = "); printr(yvec);

  for (int i = 0; i < 3; ++i)
  {
    assert(fabs(yvec[i] - expected_yvec[i]) < 1e-5 && "Assertion failed for yvec[i]");
  }
}

// BLAS-style matrix-vector multiplication with a 3D matrix (first two dimensions flattened)
__device__ void gemv312(const float *A, const float *x, float *y,
                        size_t dim1, size_t dim2, size_t dim3,
                        float alpha = 1.0f, float beta = 0.0f,
                        bool transpose = false)
{
  
  // size_t flattened_rows = dim1 * dim2;
  // size_t cols = dim3;
  // // printf("gemv312 called with A(%i %i %i) @ x(%i) + y(%i %i)\n", dim1, dim2, dim3, )
  // if (!transpose)
  // {
  //   // Flatten the first two dimensions and call gemv
  //   gemv(A, x, y, flattened_rows, cols, alpha, beta, false);
  // }
  // else
  // {
  //   // Transpose case: treat A as if it were transposed
  //   gemv(A, x, y, cols, flattened_rows, alpha, beta, true);
  // }
}

__global__ void test_gemv312()
{
  // Define a small 3D matrix A (2x2x3), vector x (3), and vector y (4)
  const int dim1 = 2, dim2 = 2, dim3 = 3;
  const int flattened_rows = dim1 * dim2;
  float A[flattened_rows * dim3] = {
      1.0f, 2.0f, 3.0f,
      4.0f, 5.0f, 6.0f,
      7.0f, 8.0f, 9.0f,
      10.0f, 11.0f, 12.0f};
  float x[dim3] = {1.0f, 1.0f, 1.0f};
  float y[flattened_rows] = {0.0f, 0.0f, 0.0f, 0.0f};

  // Call the __device__ gemv312 function directly
  gemv312(A, x, y, dim1, dim2, dim3);

  // Print the result
  printf("Result y: [%f, %f, %f, %f]\n", y[0], y[1], y[2], y[3]);

  // Expected results
  float expected_y[flattened_rows] = {6.0f, 15.0f, 24.0f, 33.0f};

  // Assert the results
  for (int i = 0; i < flattened_rows; ++i)
  {
    assert(fabs(y[i] - expected_y[i]) < 1e-5 && "Assertion failed for y[i]");
  }
}

/*
subroutine wrap_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call wrap_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine wrap_dgemv321
*/
__device__ void gemv321(const float *A, const float *x, float *y,
  size_t dim1, size_t dim2, size_t dim3,
  float alpha = 1.0f, float beta = 0.0f,
  bool transpose = false)
{
  //size_t rows = dim1;
  //size_t flattened_cols = dim2 * dim3;
//
  //if (!transpose)
  //{
  //  // Flatten the last two dimensions and call gemv
  //  gemv(A, x, y, rows, flattened_cols, alpha, beta, false);
  //}
  //else
  //{
  //  // Transpose case: treat A as if it were transposed
  //  gemv(A, x, y, flattened_cols, rows, alpha, beta, true);
  //}
}

/*
subroutine wrap_dgemv422(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)*size(amat, 4)) => amat
   xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   call wrap_gemv(aptr, xptr, yptr, alpha, beta, tra)
end subroutine wrap_dgemv422
*/

__device__ void gemv422(const float *A, const float *x, float *y,
                        size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                        float alpha = 1.0f, float beta = 0.0f,
                        bool transpose = false)
{
  // size_t flattened_rows = dim1 * dim2;
  // size_t flattened_cols = dim3 * dim4;

  // if (!transpose)
  // {
  //   // Flatten the first two dimensions and the last two dimensions, then call gemv
  //   gemv(A, x, y, flattened_rows, flattened_cols, alpha, beta, false);
  // }
  // else
  // {
  //   // Transpose case: treat A as if it were transposed
  //   gemv(A, x, y, flattened_cols, flattened_rows, alpha, beta, true);
  // }
}

__global__ void test_gemv422()
{
  // Define a small 4D matrix A (2x2x2x3), vector x (6), and vector y (4)
  const int dim1 = 2, dim2 = 2, dim3 = 2, dim4 = 3;
  const int flattened_rows = dim1 * dim2;
  const int flattened_cols = dim3 * dim4;
  float A[flattened_rows * flattened_cols] = {
      1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f,
      7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f,
      13.0f, 14.0f, 15.0f, 16.0f, 17.0f, 18.0f,
      19.0f, 20.0f, 21.0f, 22.0f, 23.0f, 24.0f};
  float x[flattened_cols] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  float y[flattened_rows] = {0.0f, 0.0f, 0.0f, 0.0f};

  // Call the __device__ gemv422 function directly
  gemv422(A, x, y, dim1, dim2, dim3, dim4);

  // Print the result
  printf("Result y: [%f, %f, %f, %f]\n", y[0], y[1], y[2], y[3]);

  // Expected results
  float expected_y[flattened_rows] = {21.0f, 57.0f, 93.0f, 129.0f};

  // Assert the results
  for (int i = 0; i < flattened_rows; ++i)
  {
    assert(fabs(y[i] - expected_y[i]) < 1e-5 && "Assertion failed for y[i]");
  }
}

void test_blas()
{
  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // test_gemv_0<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());


  cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  test_gemv_1<<<1, 1>>>();
  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());
  
  
  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // // test_gemv_1<<<1, 1>>>();

  // test_gemv312<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());

  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // // test_gemv_1<<<1, 1>>>();
  // test_gemv422<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());
  
}
