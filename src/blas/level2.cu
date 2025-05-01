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

__global__
void test_gemv312_1()
{
  printf("========================================\n");
  float amat[2][1][3] {0}; arange(amat);
  float xvec[2] {0, 1};
  float yvec[1][3] {0};
  const float alpha = 1.0f;
  const float beta = 1.0f;
  const bool trans = true;

  printf("\n amat = "); printr(amat);
  printf("\n xvec = "); printr(xvec); 
  gemv312(amat, xvec, yvec, alpha, beta, trans);
  printf("========== result =========");

  // Expected results
  float expected_yvec[1][3] {0, 1, 2};
  printf("\n result = "); printr(yvec);

  assert_isclose(expected_yvec, yvec);
}

__global__
void test_gemv312_2()
{
  printf("========================================\n");
  float amat[3][9][9]   {  {    {   0.04846839,   -0.01319809,    0.00888869,   -0.00461266,   -0.00000000,   -0.00000255,   -0.00000567,   -0.04582563,    0.10322431},    {   0.01319809,    0.25033903,    0.28724288,   -0.00061147,   -0.00015417,   -0.06118460,    0.01184802,    0.00000020,    0.00000007},    {  -0.00888869,   -0.28724288,   -0.29827035,   -0.00000000,   -0.00000000,   -0.00000000,   -0.00030407,   -0.00190949,    0.00007478},    {   0.00461266,    0.00061147,    0.00000000,   -0.00076683,   -0.00599096,   -0.00000000,    0.00000000,    0.00000000,    0.00000000},    {   0.00000000,    0.00015417,    0.00000000,    0.00599096,    0.00614513,   -0.00000000,    0.00000000,    0.00000000,    0.00000000},    {   0.00000255,    0.06118460,    0.00000000,    0.00000000,    0.00000000,    0.06118715,    0.00000000,    0.00000000,    0.00000000},    {   0.00000567,   -0.01184802,    0.00030407,   -0.00000000,   -0.00000000,   -0.00000000,   -0.01153828,    0.00000000,    0.00000000},    {   0.04582563,   -0.00000020,    0.00190949,   -0.00000000,   -0.00000000,   -0.00000000,   -0.00000000,    0.04773492,    0.00000000},    {  -0.10322431,   -0.00000007,   -0.00007478,   -0.00000000,   -0.00000000,   -0.00000000,   -0.00000000,   -0.00000000,   -0.10329916}  },  {    {  -0.00989226,    0.00530896,   -0.04324786,    0.00731743,    0.00000000,    0.00000014,   -0.00000001,   -0.04842578,    0.06915485},    {  -0.00530896,   -0.44461188,   -0.40465906,    0.03913484,    0.00112959,   -0.02218241,   -0.05272570,   -0.00000019,    0.00000001},    {   0.04324786,    0.40465906,    0.45238862,    0.00000000,    0.00000000,    0.00000000,    0.00035841,    0.00384499,    0.00027829},    {  -0.00731743,   -0.03913484,   -0.00000000,   -0.03527131,    0.01118097,   -0.00000000,   -0.00000000,   -0.00000000,   -0.00000000},    {  -0.00000000,   -0.00112959,   -0.00000000,   -0.01118097,   -0.01231056,   -0.00000000,   -0.00000000,   -0.00000000,   -0.00000000},    {  -0.00000014,    0.02218241,   -0.00000000,    0.00000000,    0.00000000,    0.02218226,   -0.00000000,   -0.00000000,    0.00000000},    {   0.00000001,    0.05272570,   -0.00035841,    0.00000000,    0.00000000,    0.00000000,    0.05236730,   -0.00000000,    0.00000000},    {   0.04842578,    0.00000019,   -0.00384499,    0.00000000,    0.00000000,    0.00000000,    0.00000000,    0.04458098,    0.00000000},    {  -0.06915485,   -0.00000001,   -0.00027829,    0.00000000,    0.00000000,   -0.00000000,   -0.00000000,   -0.00000000,   -0.06943315}  },  {    {  -0.11942402,    0.00894437,    0.01865741,    0.00367547,    0.00000000,    0.00000041,    0.00000983,   -0.15660350,    0.00589198},    {  -0.00894437,    0.06929022,   -0.01479849,    0.00432650,   -0.00040069,   -0.02813573,    0.11724347,   -0.00000044,   -0.00000002},    {  -0.01865741,    0.01479849,   -0.00891470,    0.00000000,   -0.00000000,   -0.00000000,    0.00024444,   -0.00522030,   -0.00007991},    {  -0.00367547,   -0.00432650,   -0.00000000,   -0.02915783,   -0.02115586,   -0.00000000,    0.00000000,   -0.00000000,   -0.00000000},    {  -0.00000000,    0.00040069,    0.00000000,    0.02115586,    0.02155655,    0.00000000,    0.00000000,   -0.00000000,   -0.00000000},    {  -0.00000041,    0.02813573,    0.00000000,    0.00000000,   -0.00000000,    0.02813532,    0.00000000,   -0.00000000,   -0.00000000},    {  -0.00000983,   -0.11724347,   -0.00024444,   -0.00000000,   -0.00000000,   -0.00000000,   -0.11749773,   -0.00000000,   -0.00000000},    {   0.15660350,    0.00000044,    0.00522030,    0.00000000,    0.00000000,    0.00000000,    0.00000000,    0.16182424,    0.00000000},    {  -0.00589198,    0.00000002,    0.00007991,    0.00000000,    0.00000000,    0.00000000,    0.00000000,   -0.00000000,   -0.00581205}  }};
  float xvec[9] = {   0.00000000,    0.00000000,    0.00000000,   -0.00000000,   -0.00000003,   -0.00000001,   -0.00000001,   -0.00000002,   -0.00000002};
  float yvec[3][9] = { {  -0.02694907,    0.00810465,   -0.00265063,   -0.02681478,    0.02736383,    0.05639220,   -0.00487420,    0.01200702,   -0.04257901},  {   0.02721153,   -0.02968238,    0.01284230,    0.01578207,   -0.05118861,    0.02033170,    0.02054887,    0.01268855,   -0.02853403},  {   0.00812483,   -0.01700864,   -0.00553054,   -0.10054903,    0.09667157,    0.02612885,   -0.04602112,    0.04075407,   -0.00256999}};
  const float alpha = 1.0f;
  const float beta = 1.0f;
  const bool trans = false;

  printf("\n amat = "); printr(amat);
  printf("\n xvec = "); printr(xvec, "%.6f\0"); 
  gemv312(amat, xvec, yvec, alpha, beta, trans);
  printf("========== result =========");

  // Expected results
  const float expected_yvec[3][9] = { {  -0.02694908,    0.00810465,   -0.00265063,   -0.02681478,    0.02736383,    0.05639220,   -0.00487420,    0.01200702,   -0.04257901},  {   0.02721152,   -0.02968238,    0.01284230,    0.01578207,   -0.05118861,    0.02033170,    0.02054887,    0.01268855,   -0.02853403},  {   0.00812483,   -0.01700864,   -0.00553054,   -0.10054903,    0.09667157,    0.02612885,   -0.04602112,    0.04075407,   -0.00256999}};
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
  

  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // test_gemv312_0<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());


  // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  // test_gemv312_1<<<1, 1>>>();
  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());

  cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  test_gemv312_2<<<1, 1>>>();
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