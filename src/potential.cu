#include <stdio.h>
#include "potential.h"
#include "structure.h"

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