#include "mulliken.h"

__device__
void get_mulliken_shell_charges(
  const basis_type &bas,
  const float (&smat)[MAX_NAO][MAX_NAO],
  const float (&pmat)[MAX_NSPIN][MAX_NAO][MAX_NAO],
  const float (&n0sh)[MAX_NSH],
  float (&qsh)[MAX_NSPIN][MAX_NSH]
)
{
  // Initialize qsh to zero
  for (int spin = 0; spin < MAX_NSPIN; spin++) {
    for (int ish = 0; ish < MAX_NSH; ish++) {
      qsh[spin][ish] = 0.0;
    }
  }

  // Compute Mulliken shell charges
  for (int spin = 0; spin < MAX_NSPIN; spin++) {
    for (int iao = 0; iao < bas.nao; iao++) {
      float pao = 0.0;
      for (int jao = 0; jao < bas.nao; jao++) {
        pao += pmat[spin][iao][jao] * smat[iao][jao];
      }
      qsh[spin][bas.ao2sh[iao]] -= pao;
    }
  }

  updown_to_magnet(qsh);

  for (int iao = 0; iao < bas.nao; iao++)
  {
    qsh[0][iao] += n0sh[iao];
  }
}


/*
subroutine get_mulliken_atomic_multipoles(bas, mpmat, pmat, mpat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: mpmat(:, :, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(out) :: mpat(:, :, :)

   integer :: iao, jao, spin
   real(wp) :: pao(size(mpmat, 1))

   mpat(:, :, :) = 0.0_wp
   ! $omp parallel do default(none) schedule(runtime) reduction(+:mpat) &
   ! $omp shared(bas, pmat, mpmat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao(:) = 0.0_wp
         do jao = 1, bas%nao
            pao(:) = pao + pmat(jao, iao, spin) * mpmat(:, jao, iao)
         end do
         mpat(:, bas%ao2at(iao), spin) = mpat(:, bas%ao2at(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(mpat)

end subroutine get_mulliken_atomic_multipoles
*/

/* TODO: Go FROM HERE */
/* TODO: You need to turn 3 into 6 and with templates? */
__device__
void get_mulliken_atomic_multipoles(
  const basis_type &bas,
  const float (&mpmat)[MAX_NAO][MAX_NAO][3],
  const float (&pmat)[MAX_NSPIN][MAX_NAO][MAX_NAO],
  float (&mpat)[MAX_NSPIN][MAX_NAT][3]
)
{
  // Initialize mpat to zero
  for (int spin = 0; spin < MAX_NSPIN; spin++) {
    for (int iat = 0; iat < MAX_NAT; iat++) {
      for (int dim = 0; dim < 3; dim++) {
        mpat[spin][iat][dim] = 0.0f;
      }
    }
  }

  // Compute Mulliken atomic multipoles
  for (int spin = 0; spin < MAX_NSPIN; spin++) {
    for (int iao = 0; iao < bas.nao; iao++) {
      float pao[3] = {0.0f, 0.0f, 0.0f};
      for (int jao = 0; jao < bas.nao; jao++) {
        for (int dim = 0; dim < 3; dim++) {
          pao[dim] += pmat[spin][jao][iao] * mpmat[dim][jao][iao];
        }
      }
      for (int dim = 0; dim < 3; dim++) {
        mpat[spin][bas.ao2at[iao]][dim] -= pao[dim];
      }
    }
  }

  updown_to_magnet(mpat);
}