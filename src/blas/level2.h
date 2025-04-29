#ifndef BLAS_LEVEL2_H
#define BLAS_LEVEL2_H
/* 
This function does y := alpha*A*x + beta*y,
*/
// template <int N>
__device__ __host__ 
void symv(
  int N,
  const float *amat, // 1D array representing symmetric matrix
  const float (&xvec)[],//[N],
  float (&yvec)[], //[N],
  const float beta = 0.0f,
  const float alpha = 1.0f
);

// matrix-vector multiplication
template <int N, int M>
__device__ __host__
void gemv(
  const float (&A)[N][M], const float (&x)[M], float (&y)[N], 
  float alpha = 1.0f, 
  float beta = 0.0f, 
  bool transpose = false
)
{
  for (int i = 0; i < N; ++i) {
    float temp = 0.0f;
    for (int j = 0; j < M; ++j) {
      temp += transpose ? A[j][i] * x[j] : A[i][j] * x[j];
    }
    y[i] = alpha * temp + beta * y[i];
  }
}

/*

*/
template <int N, int K, int M>
__device__  __host__
void gemv( /* gemv321 */
  const float (&A)[N][K][M], const float (&x)[K][M], float (&y)[N], 
  float alpha = 1.0f, 
  float beta = 0.0f, 
  bool transpose = false)
{
  for (int i = 0; i < N; ++i) {
    float temp = 0.0f;
    for (int j = 0; j < K; ++j) {
      for (int k = 0; k < M; ++k) {
        int idx = j * M + k; // Flatten [K][M] into [K * M]
        temp += transpose ? A[j][k][i] * x[j][k] : A[i][j][k] * x[j][k];
      }
    }
    y[i] = alpha * temp + beta * y[i];
  }
}

/*
Perform matrix-vector multiplication, with the following shapes

A(M, N, M, N) @ x(M, N) + y(M, N)
Where interpret the shapes like so:

A(M*N, M*N) @ x(M*N) + y(M*N)
*/
template <int N, int M>
__device__ __host__
void gemv(
  const float (&A)[M][N][M][N], const float (&x)[M][N], float (&y)[M][N], 
  float alpha = 1.0f, 
  float beta = 0.0f, 
  bool transpose = false)
{
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      float temp = 0.0f;
      for (int k = 0; k < M; ++k) {
        for (int l = 0; l < N; ++l) {
          int idx1 = k * N + l; // Flatten [M][N] into [M * N]
          int idx2 = j * M + i; // Transpose after flattening
          temp += transpose ? A[k][l][j][i] * x[k][l] : A[i][j][k][l] * x[k][l];
        }
      }
      y[i][j] = alpha * temp + beta * y[i][j];
    }
  }
}



void test_blas();
#endif