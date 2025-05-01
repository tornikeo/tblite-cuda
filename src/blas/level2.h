#ifndef BLAS_LEVEL2_H
#define BLAS_LEVEL2_H
#include "../utils/array.h"
#include "../limits.h"
#include "../utils/gpu.h"
#include "../utils/array.h"
#include <stdio.h>
#include <cassert>
#include <curand.h>
#include <cmath>

__device__ __host__ inline 
/* 
yvec := alpha * amat(:,:) @ xvec(:) + beta * yvec(:)

where `amat` is a symmetric matrix
*/
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
__device__ __host__ inline
/*
yvec := alpha * amat(:,:) @ xvec(:) + beta * yvec(:)

where `amat` is any matrix. Shapes are checked prior.
*/
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
__device__ __host__ inline
/*
```
if transpose:
  yvec := alpha * xvec(:*:) @ amat(:*:,:)  + beta * yvec(:)
else:
  yvec := alpha * amat(:,:*:) @ xvec(:*:) + beta * yvec(:)
```

*/
void gemv321(const float *A, const float *x, float *y,
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

__device__ __host__ inline 
void gemv422(const float *A, const float *x, float *y,
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

template <int N, int M, int L, int K, int O, int P>
__device__ __host__ inline 
void gemv312(
  const float (&amat)[N][M][L], 
  const float (&xvec)[K], 
  float (&yvec)[O][P],
  float alpha, 
  float beta,
  bool transpose
)
{
  if(transpose) // (N,M*L)^T @ K + O,P
  {
    assert(N == K);
    assert(M*L == O*P);
    for (int i = 0; i < O; i++)
    {
      for (int j = 0; j < P; j++)
      {
        float temp = 0;
        for (int k = 0; k < N; k++)
        {
          temp += amat[k][i][j] * xvec[k];
        }
        yvec[i][j] = alpha * temp + beta * yvec[i][j];
      }
    }
  } else { // N*M,L @ K + O,P
    assert(L == K);
    assert(O * P == N * M);
    for (int i = 0; i < O; i++)
    {
      for (int j = 0; j < P; j++)
      {
        float temp = 0;
        for (int k = 0; k < K; k++)
        {
          temp += amat[i][j][k] * xvec[k];
        }
        yvec[i][j] = alpha * temp + beta * yvec[i][j];
      }
    }
  }
}

void test_blas();
#endif