#ifndef LAPACK_SYGVD_H
#define LAPACK_SYGVD_H
#include "../limits.h"
#include "../utils/error.h"
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
class sygvd_solver
{
  public:
  int n = MAX_NAO; 
  int ldwork = 1 + 6*MAX_NAO + 2*(MAX_NAO*MAX_NAO);
  float liwork = 3 + 5*MAX_NAO;
  float dbmat[MAX_NAO][MAX_NAO];
  float dwork[1 + 6*MAX_NAO + 2*(MAX_NAO*MAX_NAO)]; // all shapes hardcoded
  int iwork[3 + 5*MAX_NAO];
  __device__ 
  void solve(
          float (&hmat)[MAX_NAO][MAX_NAO],
    const float (&smat)[MAX_NAO][MAX_NAO],
          float (&eval)[MAX_NAO],
          error_type &err
  );
};

void test_sygvd();

#endif