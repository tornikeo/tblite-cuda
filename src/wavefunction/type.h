#ifndef WAVEFUNCTION_TYPE_H
#define WAVEFUNCTION_TYPE_H
#include "../limits.h"
/*
   type :: wavefunction_type
      !> Electronic temperature
      real(wp) :: kt = 0.0_wp
      !> Number of electrons in this wavefunction
      real(wp) :: nocc = 0.0_wp
      !> Number of unpaired electrons in this wavefunction
      real(wp) :: nuhf = 0.0_wp
      !> Number of spin channels
      integer :: nspin = 1
      !> Index of the highest occupied molecular orbitals
      integer, allocatable :: homo(:)
      !> Number of electrons
      real(wp), allocatable :: nel(:)
      !> Reference occupation number for each atom, shape: [nat]
      real(wp), allocatable :: n0at(:)
      !> Reference occupation number for each shell, shape: [nsh]
      real(wp), allocatable :: n0sh(:)

      !> Density matrix, shape: [nao, nao, spin]
      real(wp), allocatable :: density(:, :, :)
      !> Orbital coefficients, shape: [nao, nao, spin]
      real(wp), allocatable :: coeff(:, :, :)
      !> Orbital energies, eigenvalues, shape: [nao, spin]
      real(wp), allocatable :: emo(:, :)
      !> Occupation numbers, shape: [nao, spin]
      real(wp), allocatable :: focc(:, :)

      !> Number of electrons for each atom, shape: [nat, spin]
      real(wp), allocatable :: qat(:, :)
      !> Number of electrons for each shell, shape: [nsh, spin]
      real(wp), allocatable :: qsh(:, :)

      !> Atomic dipole moments for each atom, shape: [3, nat, spin]
      real(wp), allocatable :: dpat(:, :, :)
      !> Atomic quadrupole moments for each atom, shape: [5, nat, spin]
      real(wp), allocatable :: qpat(:, :, :)
   end type wavefunction_type
*/

/*
subroutine new_wavefunction(self, nat, nsh, nao, nspin, kt)
   type(wavefunction_type), intent(out) :: self
   integer, intent(in) :: nat
   integer, intent(in) :: nsh
   integer, intent(in) :: nao
   integer, intent(in) :: nspin
   real(wp), intent(in) :: kt

   self%nspin = nspin
   self%kt = kt

   allocate(self%homo(max(2, nspin)))
   allocate(self%nel(max(2, nspin)))

   allocate(self%n0at(nat))
   allocate(self%n0sh(nsh))

   allocate(self%density(nao, nao, nspin))
   allocate(self%coeff(nao, nao, nspin))
   allocate(self%emo(nao, nspin))
   allocate(self%focc(nao, nspin))

   allocate(self%qat(nat, nspin))
   allocate(self%qsh(nsh, nspin))

   allocate(self%dpat(3, nat, nspin))
   allocate(self%qpat(6, nat, nspin))

   self%qat(:, :) = 0.0_wp
   self%qsh(:, :) = 0.0_wp
   self%dpat(:, :, :) = 0.0_wp
   self%qpat(:, :, :) = 0.0_wp
end subroutine new_wavefunction

*/

typedef struct
{
  float kt = 0;
  float nocc = 0;
  float nuhf = 0;
  int nspin = 1;
  int homo[2];
  float nel[2];
  float n0at[MAX_NAT] = {0};
  float n0sh[MAX_NSH] = {0};
  float density[MAX_NAO][MAX_NAO][MAX_NSPIN] = {0};
  float coeff[MAX_NAO][MAX_NAO][MAX_NSPIN] = {0};
  float emo[MAX_NAO][MAX_NSPIN] = {0};
  float focc[MAX_NAO][MAX_NSPIN] = {0};
  float qat[MAX_NAT][MAX_NSPIN] = {0};
  float qsh[MAX_NSH][MAX_NSPIN] = {0};
  float dpat[3][MAX_NAT][MAX_NSPIN] = {0};
  float qpat[5][MAX_NAT][MAX_NSPIN] = {0};
} wavefunction_type;

/*
!> Split an real occupation number into alpha and beta space.
!>
!> This routine does not perform any checks on the condition
!> ``mod(nocc, 2) == 0 .eqv. mod(nuhf, 2) == 0`` and will yield fractional
!> occupations in case those condtions are not fullfilled.
!> However, it will avoid creating negative occupation numbers.
subroutine get_alpha_beta_occupation(nocc, nuhf, nalp, nbet)
   real(wp), intent(in) :: nocc
   real(wp), intent(in) :: nuhf
   real(wp), intent(out) :: nalp
   real(wp), intent(out) :: nbet

   real(wp) :: ntmp, diff

   ! make sure we cannot get a negative occupation here
   diff = min(nuhf, nocc)
   ntmp = nocc - diff

   nalp = ntmp / 2 + diff
   nbet = ntmp / 2
end subroutine get_alpha_beta_occupation
*/

__device__ void get_alpha_beta_occupation(
    float nocc, float nuhf, float &nalp, float &nbet);
#endif