#include "h0.h"

/*subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Occupation number
      real(wp), intent(out) :: nocc
      !> Reference occupation for each atom
      real(wp), intent(out) :: n0at(:)
      !> Reference occupation for each shell
      real(wp), intent(out) :: n0sh(:)

      integer :: iat, ish, izp, ii

      nocc = -mol%charge
      n0at(:) = 0.0_wp
      n0sh(:) = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            nocc = nocc + h0%refocc(ish, izp)
            n0at(iat) = n0at(iat) + h0%refocc(ish, izp)
            n0sh(ii+ish) = n0sh(ii+ish) + h0%refocc(ish, izp)
         end do
      end do

   end subroutine get_occupation
   */
__device__
void get_occupation(
  const structure_type &mol,
  const basis_type &bas,
  const tb_hamiltonian &h0,
  float &nocc,
  float (&n0at)[MAX_NAT],
  float (&n0sh)[MAX_NSH]
) {
  int iat, ish, izp, ii;

  nocc = -mol.charge;
  for (int i = 0; i < MAX_NAT; ++i) {
    n0at[i] = 0.0f;
  }
  for (int i = 0; i < MAX_NSH; ++i) {
    n0sh[i] = 0.0f;
  }

  for (iat = 0; iat < mol.nat; ++iat) {
    izp = mol.id[iat];
    ii = bas.ish_at[iat];
    for (ish = 0; ish < bas.nsh_id[izp]; ++ish) {
      nocc += h0.refocc[ish][izp];
      n0at[iat] += h0.refocc[ish][izp];
      n0sh[ii + ish] += h0.refocc[ish][izp];
    }
  }
}