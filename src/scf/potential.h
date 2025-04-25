#ifndef SCF_POTENTIAL_H
#define SCF_POTENTIAL_H
#include "../limits.h"
#include "../integral/type.h"
#include "../structure.h"
#include "../basis/type.h"


// !> Container for density dependent potential-shifts
// type :: potential_type
//    !> Atom-resolved charge-dependent potential shift
//    real(wp), allocatable :: vat(:, :)
//    !> Shell-resolved charge-dependent potential shift
//    real(wp), allocatable :: vsh(:, :)
//    !> Orbital-resolved charge-dependent potential shift
//    real(wp), allocatable :: vao(:, :)

//    !> Atom-resolved dipolar potential
//    real(wp), allocatable :: vdp(:, :, :)
//    !> Atom-resolved quadrupolar potential
//    real(wp), allocatable :: vqp(:, :, :)
// contains
class potential_type {
public:
  // Atom-resolved charge-dependent potential shift
  float vat[MAX_NSPIN][MAX_NAT];
  // Shell-resolved charge-dependent potential shift
  float vsh[MAX_NSPIN][MAX_NSH];
  // Orbital-resolved charge-dependent potential shift
  float vao[MAX_NSPIN][MAX_NAO];

  // Atom-resolved dipolar potential
  float vdp[MAX_NSPIN][MAX_NAT][3];
  // Atom-resolved quadrupolar potential
  float vqp[MAX_NSPIN][MAX_NAT][6];

  // Member function to reset all values to 0.0
  __device__ void reset();
};

__device__
void new_potential(
  potential_type &self, 
  const structure_type &mol, 
  const basis_type &bas,
  const int nspin
);


/*
!> Expand an atom-resolved potential shift to a shell-resolved potential shift
subroutine add_vat_to_vsh(bas, vat, vsh)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Atom-resolved charge-dependent potential shift
   real(wp), intent(in) :: vat(:, :)
   !> Shell-resolved charge-dependent potential shift
   real(wp), intent(inout) :: vsh(:, :)

   integer :: iat, ish, ii, spin

   ! $omp parallel do schedule(runtime) collapse(2) default(none) &
   ! $omp reduction(+:vsh) shared(bas, vat) private(spin, ii, ish, iat)
   do spin = 1, size(vat, 2)
      do iat = 1, size(vat, 1)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_at(iat)
            vsh(ii+ish, spin) = vsh(ii+ish, spin) + vat(iat, spin)
         end do
      end do
   end do
end subroutine add_vat_to_vsh
*/

__device__ 
void add_pot_to_h1(
  const basis_type &bas, 
  const integral_type &ints, 
  potential_type &pot, 
  float (&h1)[MAX_NSPIN][MAX_NAO][MAX_NAO]
);

#endif