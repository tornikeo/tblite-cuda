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
  printf("Cholesky, before");
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
  printf("Cholesky, after");
  printr(n,n,b);
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
__device__
void two_sided_trsm_upper(const int n, float* a, const int lda, float* u, const int ldu) {
  // Perform A := U^{-T} * A * U^{-1}
  // Step 1: Solve Uᵀ X = A  -> overwrite A with X
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j; ++i) {
      float sum = a[i * lda + j];
      for (int k = 0; k < i; ++k) {
        sum -= u[k * ldu + i] * a[k * lda + j];
      }
      a[i * lda + j] = sum / u[i * ldu + i];
    }
  }
  // Step 2: Solve X U = A'  -> overwrite A with A'
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      float sum = a[i * lda + j];
      for (int k = 0; k < j; ++k) {
        sum -= a[i * lda + k] * u[k * ldu + j];
      }
      a[i * lda + j] = sum / u[j * ldu + j];
      if (i != j) {
        a[j * lda + i] = a[i * lda + j]; // Symmetrize
      }
    }
  }
}
  
__device__
void jacobi_eigenvalue(const int n, float* a, const int lda, float* w, float* v, const int ldv, int* info) {
  // Very naive Jacobi method for small matrices
  const int max_sweeps = 50;
  const float tol = 1e-6f;
  
  // Initialize V = I
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      v[i * ldv + j] = (i == j) ? 1.0f : 0.0f;
    }
  }

  for (int sweep = 0; sweep < max_sweeps; ++sweep) {
    bool converged = true;
    for (int p = 0; p < n - 1; ++p) {
      for (int q = p + 1; q < n; ++q) {
        float app = a[p * lda + p];
        float aqq = a[q * lda + q];
        float apq = a[p * lda + q];
        if (fabsf(apq) > tol * sqrtf(app * aqq)) {
          converged = false;
          float phi = 0.5f * atanf(2.0f * apq / (aqq - app));
          float c = cosf(phi);
          float s = sinf(phi);
          
          // Rotate rows and columns p and q
          for (int k = 0; k < n; ++k) {
            float akp = a[k * lda + p];
            float akq = a[k * lda + q];
            a[k * lda + p] = c * akp - s * akq;
            a[k * lda + q] = s * akp + c * akq;
          }
          for (int k = 0; k < n; ++k) {
            float apk = a[p * lda + k];
            float aqk = a[q * lda + k];
            a[p * lda + k] = c * apk - s * aqk;
            a[q * lda + k] = s * apk + c * aqk;
          }
          
          // Rotate eigenvectors
          for (int k = 0; k < n; ++k) {
            float vkp = v[k * ldv + p];
            float vkq = v[k * ldv + q];
            v[k * ldv + p] = c * vkp - s * vkq;
            v[k * ldv + q] = s * vkp + c * vkq;
          }
        }
      }
    }
    if (converged) break;
  }

  // Extract eigenvalues from the diagonal
  for (int i = 0; i < n; ++i) {
    w[i] = a[i * lda + i];
  }
  
  *info = 0;
}
  
/* A horrifyingly inefficient way to solve a generalized eigenvalue
  problem. TODO: priority medium, please use cuSOLVER.
*/
__device__
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
{
  // Assume itype == 1, jobz == 'v', uplo == 'u'
  if (n <= 0 || lda < n || ldb < n) {
    *info = -1;
    printf("Error: Invalid matrix dimensions or parameters (n: %d, lda: %d, ldb: %d)\n", n, lda, ldb);
    assert(false && "Bad set of parameters n lda ldb");
    return;
  }

  // Step 1: Cholesky factorization of B
  cholesky_upper(n, b, ldb, info);
  if (*info != 0) {
    return;
  }

  // Step 2: Reduce to standard eigenproblem
  two_sided_trsm_upper(n, a, lda, b, ldb);

  // Step 3: Solve standard eigenproblem (naive Jacobi method)
  jacobi_eigenvalue(n, a, lda, w, a, lda, info);
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
  ssygvd(1, 'v', 'u', n, A, lda, B, ldb, W, Work, 1, IWork, 1, &Info);

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
  ssygvd(1, 'v', 'u', n, &A[0][0], lda, &B[0][0], ldb, &lambda[0], Work, 1, IWork, 1, &Info);

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
  ssygvd(
    1, 'v', 'u', n, 
    &hmat[0][0], n, 
    &dbmat[0][0], n, 
    &eval[0], 
    &dwork[0], ldwork,
    &iwork[0], liwork, 
    &info
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
    cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
    test_cholesky_upper_0<<<1, 1>>>();
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());


    // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
    // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
    // test_ssygvd_0<<<1, 1>>>();
    // gpuErrchk(cudaPeekAtLastError());
    // gpuErrchk(cudaDeviceSynchronize());

    // cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
    // cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
    // test_ssygvd_1<<<1, 1>>>();
    // gpuErrchk(cudaPeekAtLastError());
    // gpuErrchk(cudaDeviceSynchronize());
  }
}
