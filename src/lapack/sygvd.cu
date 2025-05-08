#include "sygvd.h"
#include <stdio.h>
#include <cassert>
#include <math.h>
#include "../utils/error.h"
#include "../utils/gpu.h"
#include "../utils/array.h"

/*
pure subroutine ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssygvd
  */


__device__ 
void cholesky_upper(const int n, float* b, const int ldb, int* info) {
  // Cholesky factorization: B = Uᵀ U
  // printf("Cholesky, before");
  printr(n,n,b);
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < j; ++k) {
      b[k * ldb + j] /= b[k * ldb + k];
      for (int i = k + 1; i <= j; ++i) {
        b[i * ldb + j] -= b[k * ldb + i] * b[k * ldb + j];
      }
    }
    if (b[j * ldb + j] <= 0.0f) {
      *info = -2; // Not positive definite
      return;
    }
    b[j * ldb + j] = sqrtf(b[j * ldb + j]);
  }
  // printf("Cholesky, after");
  // printr(n,n,b);
}

__global__ 
void test_cholesky_upper_0()
{
  float b[3][3] = {
    {10.00, 2.00, 3.00},
    {2.00, 10.00, 5.00},
    {3.00, 5.00, 10.00}};
  float tol = 1e-2f;
  bool test_passed = true;
  float expected_b[3][3] = {
    {3.16227766f, 0.63245553f, 0.94868330f},
    {2.00f,      3.09838668f, 1.42f},
    {3.00f,      5.00f,      2.6615}};
  int n = 3;
  int info = 0;
  cholesky_upper(3, &b[0][0], 3, &info);
  assert_isclose(expected_b, b, 1e-2);
}

/* A horrifyingly inefficient way to solve a generalized eigenvalue
  problem. TODO: priority medium, please use cuSOLVER.
*/
__device__ 
int ssygv(int itype, char jobz, char uplo, int n, 
  float* A, int ldA, float* B, int ldB, float* w, float* work, int lwork) {
  // Step 1: Validate input
  // if (n <= 0 || !A || !B || !w || ldA < n || ldB < n || lwork < n) {
  //     return -1; // Invalid input
  // }

  // // Step 2: Perform Cholesky factorization of B (upper triangular)
  // int info = cholesky_upper(n, B, ldB);
  // if (info != 0) {
  //     return info; // Cholesky factorization failed
  // }

  // // Step 3: Transform A into a standard eigenvalue problem
  // transform_matrix(n, A, ldA, B, ldB);

  // // Step 4: Compute eigenvalues and eigenvectors of the transformed matrix
  // info = symmetric_eigen_decomp(n, A, ldA, w, work, lwork);
  // if (info != 0) {
  //     return info; // Eigenvalue decomposition failed
  // }

  // // Step 5: Backtransform eigenvectors to original problem
  // backtransform_eigenvectors(n, B, ldB, work, n);

  // // Step 6: Copy back eigenvectors to A (if needed)
  // for (int i = 0; i < n; ++i) {
  //     for (int j = 0; j < n; ++j) {
  //         A[i * ldA + j] = work[i * n + j];
  //     }
  // }

  return 0; // Success
}

__global__
void test_ssygvd_0()
{
  printf("========================================\n");
  printf("Test in %s at %s:%d\n", __func__, __FILE__, __LINE__);
  printf("========================================\n");
  // Example small 2x2 matrices (easy to verify manually)
  const int n = 2;
  const int lda = n;
  const int ldb = n;
  float A[lda * n] = { 4.0f, 1.0f,
                        1.0f, 3.0f };
  float B[ldb * n] = { 2.0f, 0.5f,
                        0.5f, 1.5f };
  float W[n]; // eigenvalues
  float Work[1]; // dummy
  int IWork[1];  // dummy
  int Info = 0;

  // Call the naive ssygvd
//   __device__ 
// int ssygv(char jobz, char uplo, int n, 
//   float* A, int ldA, float* B, int ldB, float* w, float* work, int lwork)
  ssygv(1, 'v', 'u', n, A, lda, B, ldb, W, Work, 1);

  if (Info != 0) {
    printf("ssygvd failed with info = %d\n", Info);
    return;
  }

  printf("========== result =========");
  
  // Output eigenvalues
  printf("Eigenvalues:\n");
  for (int i = 0; i < n; ++i) {
    printf("W[%d] = %f\n", i, W[i]);
  }
  
  // Output eigenvectors (stored in A)
  printf("Eigenvectors (columns):\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%f ", A[i * lda + j]);
    }
    printf("\n");
  }
  
  // assert_isclose(expected_yvec, yvec);
}

__global__
void test_ssygvd_1()
{
  printf("========================================\n");
  printf("Test in %s at %s:%d\n", __func__, __FILE__, __LINE__);
  printf("========================================\n");
  // Example small 2x2 matrices (easy to verify manually)
  const int n = 3;
  const int lda = n;
  const int ldb = n;
  /*
    *       | 3.5 0.5 0 |
    *   A = | 0.5 3.5 0 |
    *       | 0   0   2 |
    *
    *       | 10  2   3 |
    *   B = | 2  10   5 |
    *       | 3   5  10 |
    */
  float A[3][3] = {3.5, 0.5, 0, 0.5, 3.5, 0.0, 0.0, 0.0, 2.0};
  float B[3][3] = {10.0, 2.0, 3.0, 2.0, 10.0, 5.0, 3.0, 5.0, 10.0};
  float expected_lambda[] = {0.158660256604, 0.370751508101882, 0.6}; // Known eigenvalues
  float lambda[n] = {0};
  float Work[1]; // dummy
  int IWork[1];  // dummy
  int Info = 0;

  printf("A = \n"); printr(A); printf("\n");
  printf("B = \n"); printr(B); printf("\n");
  // Call the naive ssygvd
  
  // ssygv(1, 'v', 'u', n, &A[0][0], lda, &B[0][0], ldb, &lambda[0], Work, 1, IWork, 1, &Info);

  ssygv(1, 'v', 'u', n, &A[0][0], lda, &B[0][0], ldb, &lambda[0], Work, 1);

  if (Info != 0) {
    printf("ssygvd failed with info = %d\n", Info);
    return;
  }

  printf("========== result =========\n");
  
  // Output eigenvalues
  printf("Eigenvalues:\n");
  for (int i = 0; i < n; ++i) {
    printf("W[%d] = %f\n", i, lambda[i]);
  }

  // Output eigenvectors (stored in A)
  printf("Eigenvectors (columns):\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%f ", A[i * lda + j]);
    }
    printf("\n");
  }

  assert_isclose(expected_lambda, lambda);
}

__device__
void sygvd_solver::solve(
        float (&hmat)[MAX_NAO][MAX_NAO],
  const float (&smat)[MAX_NAO][MAX_NAO],
        float (&eval)[MAX_NAO],
        error_type &err /* TODO: Priority medium. Err is not used. */
)
{
  n = MAX_NAO;
  /* n = MAX_NAO; */
  /* 
  if (self%n == 0) then
      self%n = size(hmat, 1)
   end if
   if (.not.allocated(self%dwork)) then
      allocate(self%dwork(1 + 6*self%n + 2*self%n**2))
   end if
   if (.not.allocated(self%iwork)) then
      allocate(self%iwork(3 + 5*self%n))
   end if
  */
  const int ldwork = sizeof(dwork)/sizeof(dwork[0]);
  const int liwork = sizeof(iwork)/sizeof(iwork[0]);
  /*
  void ssygvd(
  const int itype,
  const char jobz,
  const char uplo,
  const int n,
  float* a,
  const int lda,
  float* b,
  const int ldb,
  float* w,
  float* work, // Not used here
  const int lwork,
  int* iwork,  // Not used here
  const int liwork,
  int* info
)
  */
  for (size_t i = 0; i < MAX_NAO; i++)
  {
    for (size_t j = 0; j < MAX_NAO; j++)
    {
      dbmat[i][j] = smat[i][j];
    }
  }

  int info = 0;
  /* ssygvd(1, 'v', 'u', n, A, lda, B, ldb, W, Work, 1, IWork, 1, &Info); */
  ssygv(
    1, 'v', 'u', n, 
    &hmat[0][0], n, 
    &dbmat[0][0], n, 
    &eval[0], 
    &dwork[0], ldwork
  );


  printf("Solver called.\n");
  if (info != 0) {
    printf("sygvd_solver failed with info = %d\n", info);
    switch (info) {
      case -1:
      printf("Error: Invalid input parameters.\n");
      break;
      case -2:
      printf("Error: Matrix B is not positive definite.\n");
      break;
      default:
      printf("Error: Unknown failure with info = %d.\n", info);
      break;
    }
    return;
  }

  printf("Eigenvalues:\n");
  for (int i = 0; i < n; ++i) {
    printf("eval[%d] = %f\n", i, eval[i]);
  }

  printf("Eigenvectors (columns):\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%f ", hmat[i][j]);
    }
    printf("\n");
  }
}

void test_sygvd()
{
  {
    // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
    // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
    // test_cholesky_upper_0<<<1, 1>>>();
    // gpuErrchk(cudaPeekAtLastError());
    // gpuErrchk(cudaDeviceSynchronize());


    cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
    test_ssygvd_0<<<1, 1>>>();
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
    // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
    // test_ssygvd_1<<<1, 1>>>();
    // gpuErrchk(cudaPeekAtLastError());
    // gpuErrchk(cudaDeviceSynchronize());
  }
}
