#include "type.h"

__device__ void get_alpha_beta_occupation(
  float nocc, float nuhf, float &nalp, float &nbet)
{
  float ntmp, diff;

  // Ensure we cannot get a negative occupation here
  diff = fminf(nuhf, nocc);
  ntmp = nocc - diff;

  nalp = ntmp / 2.0f + diff;
  nbet = ntmp / 2.0f;
}

/*
subroutine get_density_matrix(focc, coeff, pmat)
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(out) :: pmat(:, :)

   real(wp), allocatable :: scratch(:, :)
   integer :: iao, jao

   allocate(scratch(size(pmat, 1), size(pmat, 2)))
   ! $omp parallel do collapse(2) default(none) schedule(runtime) &
   ! $omp shared(scratch, coeff, focc, pmat) private(iao, jao)
   do iao = 1, size(pmat, 1)
      do jao = 1, size(pmat, 2)
         scratch(jao, iao) = coeff(jao, iao) * focc(iao)
      end do
   end do
   call gemm(scratch, coeff, pmat, transb='t')
end subroutine get_density_matrix
*/

__device__
void get_density_matrix(
  const float (&focc)[MAX_NAO],
  const float (&coeff)[MAX_NAO][MAX_NAO],
  float (&pmat)[MAX_NAO][MAX_NAO]
)
{
  for (int iao = 0; iao < MAX_NAO; ++iao) {
    for (int jao = 0; jao < MAX_NAO; ++jao) {
      float scratch = 0.0f;
      for (int k = 0; k < MAX_NAO; ++k) {
        scratch += coeff[iao][k] * focc[k];
      }
      pmat[iao][jao] = 0.0f;
      for (int k = 0; k < MAX_NAO; ++k) {
        pmat[iao][jao] += scratch * coeff[jao][k];
      }
    }
  }
}