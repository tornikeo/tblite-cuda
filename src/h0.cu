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

__device__
void get_selfenergy(
  const tb_hamiltonian &h0,
  const int (&id)[MAX_NAT],
  const int (&ish_at)[MAX_NAT],
  const int (&nshell)[MAX_NELEM],
  const float (&cn)[MAX_NAT],
  // const float (&qat)[MAX_NAT],
  float selfenergy[MAX_NSH], // static size
  float dsedcn[MAX_NSH]     // static size
  // float dsedq[MAX_NSH],      // static size
)
{
  /* subroutine get_selfenergy(h0, id, ish_at, nshell, cn, qat, selfenergy, dsedcn, dsedq)
      type(tb_hamiltonian), intent(in) :: h0
      integer, intent(in) :: id(:)
      integer, intent(in) :: ish_at(:)
      integer, intent(in) :: nshell(:)
      real(wp), intent(in), optional :: cn(:)
      real(wp), intent(in), optional :: qat(:)
      real(wp), intent(out) :: selfenergy(:)
      real(wp), intent(out), optional :: dsedcn(:)
      real(wp), intent(out), optional :: dsedq(:)

      integer :: iat, izp, ish, ii

      selfenergy(:) = 0.0_wp
      if (present(dsedcn)) dsedcn(:) = 0.0_wp
      if (present(dsedq)) dsedq(:) = 0.0_wp
      do iat = 1, size(id)
         izp = id(iat)
         ii = ish_at(iat)
         do ish = 1, nshell(izp)
            selfenergy(ii+ish) = h0%selfenergy(ish, izp)
         end do
      end do
      if (present(cn)) then
         if (present(dsedcn)) then
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
                  dsedcn(ii+ish) = -h0%kcn(ish, izp)
               end do
            end do
         else */
  for (int i = 0; i < MAX_NSH; ++i) {
    selfenergy[i] = 0.0f;
    dsedcn[i] = 0.0f;
  }

  for (int iat = 0; iat < MAX_NAT; ++iat) {
    int izp = id[iat];
    int ii = ish_at[iat];
    for (int ish = 0; ish < nshell[izp]; ++ish) {
      selfenergy[ii + ish] = h0.selfenergy[ish][izp];
    }
  }

  for (int iat = 0; iat < MAX_NAT; ++iat) {
    int izp = id[iat];
    int ii = ish_at[iat];
    for (int ish = 0; ish < nshell[izp]; ++ish) {
      selfenergy[ii + ish] -= h0.kcn[ish][izp] * cn[iat];
      dsedcn[ii + ish] = -h0.kcn[ish][izp];
    }
  }
  /* TODO: Low priority */
  /* Other options qat, dsedq are not implemented */
}