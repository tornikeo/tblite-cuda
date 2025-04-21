#ifndef XTB_H
#define XTB_H
#include "limits.h"
#include "structure.h"
#include "basis/type.h"
#include "potential.h"
#include "coulomb.h"
#include "h0.h"
#include "wavefunction/type.h"

/*
   type :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq1(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq2(:, :)
      !> Enhancement factor to scale the Hamiltonian elements
      real(wp), allocatable :: hscale(:, :, :, :)
      !> Polynomial coefficients for distance dependent enhancement factor
      real(wp), allocatable :: shpoly(:, :)
      !> Atomic radius for polynomial enhancement
      real(wp), allocatable :: rad(:)
      !> Reference occupation numbers
      real(wp), allocatable :: refocc(:, :)
   end type tb_hamiltonian

   subroutine new_hamiltonian(self, mol, bas, spec)
      type(tb_hamiltonian), intent(out) :: self
      type(structure_type), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      class(tb_h0spec), intent(in) :: spec

      integer :: mshell

      mshell = maxval(bas%nsh_id)
      allocate(self%selfenergy(mshell, mol%nid), self%kcn(mshell, mol%nid), &
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
      call spec%get_selfenergy(mol, bas, self%selfenergy)
      call spec%get_cnshift(mol, bas, self%kcn)
      call spec%get_q1shift(mol, bas, self%kq1)
      call spec%get_q2shift(mol, bas, self%kq2)

      allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
      call spec%get_hscale(mol, bas, self%hscale)

      allocate(self%shpoly(mshell, mol%nid), self%rad(mol%nid))
      call spec%get_rad(mol, bas, self%rad)
      call spec%get_shpoly(mol, bas, self%shpoly)

      allocate(self%refocc(mshell, mol%nid))
      call spec%get_reference_occ(mol, bas, self%refocc)
   end subroutine new_hamiltonian
*/
// typedef struct {
//   float *selfenergy; // Flattened 2D array, [mshell, mol%nid]
//   float *kcn;        // Flattened 2D array, [mshell, mol%nid]
//   float *kq1;        // Flattened 2D array, [mshell, mol%nid]
//   float *kq2;        // Flattened 2D array, [mshell, mol%nid]
//   float *hscale;     // Flattened 4D array, [mshell, mshell, mol%nid, mol%nid]
//   float *shpoly;     // Flattened 2D array, [mshell, mol%nid]
//   float *rad;        // 1D array, [mol%nid]
//   float *refocc;     // Flattened 2D array, [mshell, mol%nid]
// } tb_hamiltonian;

/*   !> Container to evaluate classical repulsion interactions for the xTB Hamiltonian
   type, extends(repulsion_type) :: tb_repulsion
      !> Exponent for the repulsion interaction
      real(wp), allocatable :: alpha(:, :)
      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:, :)
      !> Scaling of the repulsion exponents
      real(wp), allocatable :: kexp(:, :)
      !> Exponent of the repulsion polynomial
      real(wp), allocatable :: rexp(:, :)
      !> Real-space cutoff
      real(wp) :: cutoff = 25.0_wp
   contains
      procedure :: get_engrad
   end type tb_repulsion*/
typedef struct {
  float alpha[MAX_NELEM][MAX_NELEM];  // Flattened 2D array, [nelem, nelem]
  float zeff[MAX_NELEM][MAX_NELEM];   // Flattened 2D array, [nelem, nelem]
  float kexp[MAX_NELEM][MAX_NELEM];   // Flattened 2D array, [nelem, nelem]
  float rexp[MAX_NELEM][MAX_NELEM];   // Flattened 2D array, [nelem, nelem]
  float cutoff = 25.0f; // Real-space cutoff
} tb_repulsion;

typedef struct 
{
/* data */

} container_cache;


// typedef struct 
// {
// /* data */
// } dispersion_type;


typedef struct
{
  // type(basis_type) :: bas
  // type(tb_hamiltonian) :: h0
  // class(ncoord_type), allocatable :: ncoord
  // type(tb_repulsion), allocatable :: repulsion
  // type(tb_coulomb), allocatable :: coulomb
  // type(halogen_correction), allocatable :: halogen
  // class(dispersion_type), allocatable :: dispersion
  // real(wp) :: mixer_damping = mixer_damping_default
  // integer :: max_iter = max_iter_default
  // logical :: save_integrals = .false.
  // !> List of additional interaction containers
  // type(container_list), allocatable :: interactions
  basis_type bas;
  tb_hamiltonian h0;
  gfn_ncoord_type ncoord;
  tb_repulsion repulsion;
  tb_coulomb coulomb; 
  // dispersion_type dispersion;
} xtb_calculator;

void xtb_test();

#endif