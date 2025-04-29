#ifndef WAVEFUNCTION_MULLIKEN_H
#define WAVEFUNCTION_MULLIKEN_H
#include "../basis/type.h"
#include "../limits.h"
#include "spin.h"
/*

subroutine get_mulliken_shell_charges(bas, smat, pmat, n0sh, qsh)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: smat(:, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(in) :: n0sh(:)
   real(wp), intent(out) :: qsh(:, :)

   integer :: iao, jao, spin
   real(wp) :: pao

   qsh(:, :) = 0.0_wp
   ! $omp parallel do default(none) collapse(2) schedule(runtime) reduction(+:qsh) &
   ! $omp shared(bas, pmat, smat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao = 0.0_wp
         do jao = 1, bas%nao
            pao = pao + pmat(jao, iao, spin) * smat(jao, iao)
         end do
         qsh(bas%ao2sh(iao), spin) = qsh(bas%ao2sh(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(qsh)
   qsh(:, 1) = qsh(:, 1) + n0sh

end subroutine get_mulliken_shell_charges

*/
__device__
void get_mulliken_shell_charges(
  const basis_type &bas,
  const float (&smat)[MAX_NAO][MAX_NAO],
  const float (&pmat)[MAX_NSPIN][MAX_NAO][MAX_NAO],
  const float (&n0sh)[MAX_NSH],
  float (&qsh)[MAX_NSPIN][MAX_NSH]
);

// template <int Dim>
// __device__
// void get_mulliken_atomic_multipoles(
//   const basis_type &bas,
//   const float (&mpmat)[MAX_NAO][MAX_NAO][Dim],
//   const float (&pmat)[MAX_NSPIN][MAX_NAO][MAX_NAO],
//   float (&mpat)[MAX_NSPIN][MAX_NAT][Dim]
// );

// __device__
// void get_mulliken_atomic_multipoles(
//   const basis_type &bas,
//   const float (&mpmat)[MAX_NAO][MAX_NAO][3],
//   const float (&pmat)[MAX_NSPIN][MAX_NAO][MAX_NAO],
//   float (&mpat)[MAX_NSPIN][MAX_NAT][3]
// );

template <int Dim> /* Template functions need to be in the header to be used*/
__device__         /* in other files */
void get_mulliken_atomic_multipoles(
  const basis_type &bas,
  const float (&mpmat)[MAX_NAO][MAX_NAO][Dim],
  const float (&pmat)[MAX_NSPIN][MAX_NAO][MAX_NAO],
  float (&mpat)[MAX_NSPIN][MAX_NAT][Dim]
)
{
  // Initialize mpat to zero
  for (int spin = 0; spin < MAX_NSPIN; spin++) {
    for (int iat = 0; iat < MAX_NAT; iat++) {
      for (int dim = 0; dim < Dim; dim++) {
        mpat[spin][iat][dim] = 0.0f;
      }
    }
  }

  // Compute Mulliken atomic multipoles
  for (int spin = 0; spin < MAX_NSPIN; spin++) {
    for (int iao = 0; iao < bas.nao; iao++) {
      float pao[6] = {0.0};
      for (int jao = 0; jao < bas.nao; jao++) {
        for (int dim = 0; dim < Dim; dim++) {
          pao[dim] += pmat[spin][jao][iao] * mpmat[dim][jao][iao];
        }
      }
      for (int dim = 0; dim < Dim; dim++) {
        mpat[spin][bas.ao2at[iao]][dim] -= pao[dim];
      }
    }
  }

  updown_to_magnet(mpat);
}

#endif