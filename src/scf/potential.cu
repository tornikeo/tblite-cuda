#include <stdio.h>
#include "potential.h"
#include "structure.h"

__device__
void potential_type::reset() {
  for (int i = 0; i < MAX_NSPIN; ++i) {
    for (int j = 0; j < MAX_NAT; ++j) {
      vat[i][j] = 0.0f;
    }
    for (int j = 0; j < MAX_NSH; ++j) {
      vsh[i][j] = 0.0f;
    }
    for (int j = 0; j < MAX_NAO; ++j) {
      vao[i][j] = 0.0f;
    }
  }

  for (int i = 0; i < MAX_NAT; ++i) {
    for (int j = 0; j < MAX_NSPIN; ++j) {
      for (int k = 0; k < 3; ++k) {
        vdp[i][j][k] = 0.0f;
      }
      for (int k = 0; k < 6; ++k) {
        vqp[i][j][k] = 0.0f;
      }
    }
  }
}

/*subroutine new_potential(self, mol, bas, nspin)
   !> Instance of the density dependent potential
   type(potential_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Number of spin channels
   integer, intent(in) :: nspin

   allocate(self%vat(mol%nat, nspin))
   allocate(self%vsh(bas%nsh, nspin))
   allocate(self%vao(bas%nao, nspin))

   allocate(self%vdp(3, mol%nat, nspin))
   allocate(self%vqp(6, mol%nat, nspin))
end subroutine new_potential
*/

__device__
void new_potential(
  potential_type &self, 
  const structure_type &mol, 
  const basis_type &bas,
  const int nspin
)
{
  // printf("Here we are!\n");
  // do nothing since arrays are stack allocated at max already
}



template <int NSPIN, int NAT, int NSH>
__device__
void add_vat_to_vsh(
  const basis_type &bas,
  const float (&vat)[NSPIN][NAT],
        float (&vsh)[NSPIN][NSH]
)
{
  for (int spin = 0; spin < NSPIN; ++spin) {
    for (int iat = 0; iat < NAT; ++iat) {
      int ii = bas.ish_at[iat];
      for (int ish = 0; ish < bas.nsh_at[iat]; ++ish) {
        vsh[spin][ii + ish] += vat[spin][iat];
      }
    }
  }
}

/*
subroutine add_vsh_to_vao(bas, vsh, vao)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(in) :: vsh(:, :)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vao(:, :)

   integer :: ish, iao, ii, spin

   ! $omp parallel do schedule(runtime) collapse(2) default(none) &
   ! $omp reduction(+:vao) shared(bas, vsh) private(ii, iao, ish)
   do spin = 1, size(vsh, 2)
      do ish = 1, size(vsh, 1)
         ii = bas%iao_sh(ish)
         do iao = 1, bas%nao_sh(ish)
            vao(ii+iao, spin) = vao(ii+iao, spin) + vsh(ish, spin)
         end do
      end do
   end do
end subroutine add_vsh_to_vao
*/

template <int NSPIN, int NSH, int NAO>
__device__
void add_vsh_to_vao(
  const basis_type &bas,
  const float (&vat)[NSPIN][NSH],
        float (&vsh)[NSPIN][NAO]
)
{
  for (int spin = 0; spin < NSPIN; ++spin) {
    for (int ish = 0; ish < NSH; ++ish) {
      int ii = bas.iao_sh[ish];
      for (int iao = 0; iao < bas.nao_sh[ish]; ++iao) {
        vsh[spin][ii + iao] += vat[spin][ish];
      }
    }
  }
}

/*
subroutine add_vao_to_h1(bas, sint, vao, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap integrals
   real(wp), intent(in) :: sint(:, :)
   !> Orbital-resolved charge-dependent potential shift
   real(wp), intent(in) :: vao(:, :)
   !> Effective Hamiltonian
   real(wp), intent(inout) :: h1(:, :, :)

   integer :: iao, jao, spin

   ! $omp parallel do collapse(3) schedule(runtime) default(none) &
   ! $omp shared(h1, bas, sint, vao) private(spin, iao, jao)
   do spin = 1, size(h1, 3)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            h1(jao, iao, spin) = h1(jao, iao, spin) &
               & - sint(jao, iao) * 0.5_wp * (vao(jao, spin) + vao(iao, spin))
         end do
      end do
   end do
end subroutine add_vao_to_h1
*/

template <int NSPIN, int NAO>
__device__
void add_vao_to_h1(
  const basis_type &bas,
  const float (&sint)[NAO][NAO],
  const float (&vao)[NSPIN][NAO],
        float (&h1)[NSPIN][NAO][NAO]
)
{
  for (int spin = 0; spin < NSPIN; ++spin) {
    for (int iao = 0; iao < NAO; ++iao) {
      for (int jao = 0; jao < NAO; ++jao) {
        h1[spin][jao][iao] -= sint[jao][iao] * 0.5f * (vao[spin][jao] + vao[spin][iao]);
      }
    }
  }
}

/*
subroutine add_vmp_to_h1(bas, mpint, vmp, h1)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Multipole integrals, multipole operator is always centered on last index
   real(wp), intent(in) :: mpint(:, :, :)
   !> Multipole potential
   real(wp), intent(in) :: vmp(:, :, :)
   !> Effective Hamiltonian
   real(wp), intent(inout) :: h1(:, :, :)

   integer :: iao, jao, iat, jat, nmp, spin

   nmp = min(size(mpint, 1), size(vmp, 1))

   ! $omp parallel do collapse(3) schedule(runtime) default(none) &
   ! $omp shared(h1, bas, mpint, vmp, nmp) private(spin, iao, jao, iat, jat)
   do spin = 1, size(h1, 3)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            iat = bas%ao2at(iao)
            jat = bas%ao2at(jao)
            h1(jao, iao, spin) = h1(jao, iao, spin) &
               & - 0.5_wp * dot_product(mpint(:nmp, jao, iao), vmp(:nmp, iat, spin)) &
               & - 0.5_wp * dot_product(mpint(:nmp, iao, jao), vmp(:nmp, jat, spin))
         end do
      end do
   end do
end subroutine add_vmp_to_h1
*/
template <int D, int NAO, int NAT, int NSPIN>
__device__
void add_vmp_to_h1( /* TODO: Is this correct? */
  const basis_type &bas,
  const float (&mpint)[NAO][NAO][D],
  const float (&vmp)[NSPIN][NAT][D],
        float (&h1)[NSPIN][NAO][NAO]
)
{
  for (int spin = 0; spin < NSPIN; ++spin) {
    for (int iao = 0; iao < NAO; ++iao) {
      for (int jao = 0; jao < NAO; ++jao) {
        int iat = bas.ao2at[iao];
        int jat = bas.ao2at[jao];
        float dot1 = 0.0f;
        float dot2 = 0.0f;
        for (int d = 0; d < D; ++d) {
          dot1 += mpint[jao][iao][d] * vmp[spin][iat][d];
          dot2 += mpint[iao][jao][d] * vmp[spin][jat][d];
        }
        h1[spin][jao][iao] -= 0.5f * (dot1 + dot2);
      }
    }
  }
}

__device__ 
void add_pot_to_h1(
  const basis_type &bas, 
  const integral_type &ints, 
  potential_type &pot, 
  float (&h1)[MAX_NSPIN][MAX_NAO][MAX_NAO]
)
{
  /*
   h1(:, :, 1) = ints%hamiltonian
   if (size(h1, 3) > 1) h1(:, :, 2:) = 0.0_wp

   call add_vat_to_vsh(bas, pot%vat, pot%vsh)
   call add_vsh_to_vao(bas, pot%vsh, pot%vao)
   call add_vao_to_h1(bas, ints%overlap, pot%vao, h1)
   call add_vmp_to_h1(bas, ints%dipole, pot%vdp, h1)
   call add_vmp_to_h1(bas, ints%quadrupole, pot%vqp, h1)

   call magnet_to_updown(h1)
  */
  // Initialize h1 with the Hamiltonian integrals

  /* h1(:, :, 1) = ints%hamiltonian */
  for (int spin = 0; spin < MAX_NSPIN; ++spin) {
    for (int i = 0; i < MAX_NAO; ++i) {
      for (int j = 0; j < MAX_NAO; ++j) {
        h1[spin][i][j] = ints.hamiltonian[i][j];
      }
    }
  }

  /* if (size(h1, 3) > 1) h1(:, :, 2:) = 0.0_wp */
  if (MAX_NSPIN > 1)
  {
    for (int spin = 1; spin < MAX_NSPIN; ++spin) {
      for (int i = 0; i < MAX_NAO; ++i) {
        for (int j = 0; j < MAX_NAO; ++j) {
          h1[spin][i][j] = 0.0f;
        }
      }
    }
  }

  add_vat_to_vsh(bas, pot.vat, pot.vsh);
  add_vsh_to_vao(bas, pot.vsh, pot.vao);
  add_vao_to_h1(bas, ints.overlap, pot.vao, h1);
  add_vmp_to_h1(bas, ints.dipole, pot.vdp, h1);
  add_vmp_to_h1(bas, ints.quadrupole, pot.vqp, h1);
}
