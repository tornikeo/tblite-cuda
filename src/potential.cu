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
  // do nothing 
}