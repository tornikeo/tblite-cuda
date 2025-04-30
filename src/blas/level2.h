#ifndef BLAS_LEVEL2_H
#define BLAS_LEVEL2_H
/* 
This function does y := alpha*avec*xvec + beta*yvec
*/

__device__
void symv(
  const float (&amat)[MAX_NSH][MAX_NSH],
  const float (&xvec)[MAX_NSH],
  float (&yvec)[MAX_NSH],
  const float alpha,
  const float beta
);

// matrix-vector multiplication
__device__  __host__
void gemv( 
  const float (&amat)[MAX_NAO][MAX_NAT],
  const float (&xvec)[MAX_NAT],
  float (&yvec)[MAX_NAO],
  const float alpha,
  const float beta,
  const bool trans
);

// matrix-vector multiplication, 3D @ 1D + 2D. e.g:
// y(3,9) := A(3,9,9) @ x(9) + y(3,9)
// y(3*9) := A(3*9,9) @ x(9) + y(3*9)
// y(27) := A(27,9) @ x(9) + y(27)
__device__ void gemv312(const float* A, const float* x, float* y, 
  size_t dim1, size_t dim2, size_t dim3, 
  float alpha = 1.0f, float beta = 0.0f, 
  bool transpose = false);

// matrix-vector multiplication, 3D @ 2D + 1D. e.g:
// y(9) := A(9,3,9) @ x(3,9) + y(9)
// y(9) := A(9,27) @ x(27) + y(9) 
__device__ void gemv321(const float *A, const float *x, float *y,
  size_t dim1, size_t dim2, size_t dim3,
  float alpha = 1.0f, float beta = 0.0f,
  bool transpose = false);

__device__ void gemv422(const float* A, const float* x, float* y, 
    size_t dim1, size_t dim2, size_t dim3, size_t dim4, 
    float alpha = 1.0f, float beta = 0.0f, 
    bool transpose = false);

void test_blas();
#endif