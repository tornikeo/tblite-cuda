#include "level2.h"


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

// __global__ 
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

void test_gemv422()
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

__global__
void test_gemv312_0()
{
  printf("========================================\n");
  float amat[2][3][2] {0}; arange(amat);
  float xvec[2] {1, 1};
  float yvec[2][3] {0};
  const float alpha = 1.0f;
  const float beta = 1.0f;
  const bool trans = false;

  printf("\n amat = "); printr(amat);
  printf("\n xvec = "); printr(xvec); 
  gemv312(amat, xvec, yvec, alpha, beta, trans);
  printf("========== result =========");

  // Expected results
  float expected_yvec[2][3] {1, 5, 9, 13, 17, 21};
  printf("\n result = "); printr(yvec);

  assert_isclose(expected_yvec, yvec);
}

void test_blas()
{
  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // test_gemv_0<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());


  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // test_gemv_1<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());
  

  cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  test_gemv312_0<<<1, 1>>>();
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