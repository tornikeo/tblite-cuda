#ifndef INTEGRAL_TYPE_H
#define INTEGRAL_TYPE_H
#include "../limits.h"

/*    public :: integral_type, new_integral

   !> Integral container to store all overlap related integrals
   type :: integral_type
      !> Effective one-electron Hamiltonian
      real(wp), allocatable :: hamiltonian(:, :)
      !> Overlap integrals
      real(wp), allocatable :: overlap(:, :)
      !> Dipole moment integrals, moment operator is centered on last index
      real(wp), allocatable :: dipole(:, :, :)
      !> Quadrupole moment integrals, moment operator is centered on last index
      real(wp), allocatable :: quadrupole(:, :, :)
   end type integral_type

contains

*/

typedef struct 
{
  float hamiltonian[MAX_NAO][MAX_NAO];
  float overlap[MAX_NAO][MAX_NAO];
  float dipole[MAX_NAO][MAX_NAO][3];
  float quadrupole[MAX_NAO][MAX_NAO][6]; /* This thing is HEAVY */
} integral_type;


/*
!> Create and allocate a new integral container storage
subroutine new_integral(self, nao)
   !> Instance of the integral container
   type(integral_type), intent(out) :: self
   !> Dimension of the integrals
   integer, intent(in) :: nao

   allocate(self%hamiltonian(nao, nao))
   allocate(self%overlap(nao, nao))
   allocate(self%dipole(3, nao, nao))
   allocate(self%quadrupole(6, nao, nao))
end subroutine new_integral
*/
__device__
void new_integral(integral_type &self, const int &nao);

#endif