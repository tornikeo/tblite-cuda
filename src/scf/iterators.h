#ifndef SCF_ITERATOR_H
#define SCF_ITERATOR_H

#include "info.h"
#include "../basis/type.h"
#include "../limits.h"
#include "../structure.h"
#include "../integral/type.h"
#include "broyden.h"
#include "../lapack/sygvd.h"
#include "../wavefunction/type.h"
#include "../coulomb.h"
#include "potential.h"

/*
function get_mixer_dimension(mol, bas, info) result(ndim)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas
   type(scf_info), intent(in) :: info
   integer :: ndim

   ndim = 0

   select case(info%charge)
   case(atom_resolved)
      ndim = ndim + mol%nat
   case(shell_resolved)
      ndim = ndim + bas%nsh
   end select

   select case(info%dipole)
   case(atom_resolved)
      ndim = ndim + 3*mol%nat
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      ndim = ndim + 6*mol%nat
   end select
end function get_mixer_dimension 
*/

__device__
// constexpr /* TODO: priority low. Figure out why constexpr doesn't work for __device__ funcs */
int get_mixer_dimension(
  const structure_type &mol, 
  const basis_type &bas,
  const scf_info &info
);


/*
subroutine next_scf(iscf, mol, bas, wfn, solver, mixer, info, coulomb, dispersion, &
      & interactions, ints, pot, cache, dcache, icache, &
      & energies, error)
   !> Current iteration count
   integer, intent(inout) :: iscf
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Convergence accelerator
   type(broyden_mixer), intent(inout) :: mixer
   !> Information on wavefunction data used to construct Hamiltonian
   type(scf_info), intent(in) :: info
   !> Container for coulombic interactions
   type(tb_coulomb), intent(in), optional :: coulomb
   !> Container for dispersion interactions
   class(dispersion_type), intent(in), optional :: dispersion
   !> Container for general interactions
   type(container_list), intent(in), optional :: interactions

   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Restart data for coulombic interactions
   type(container_cache), intent(inout) :: cache
   !> Restart data for dispersion interactions
   type(container_cache), intent(inout) :: dcache
   !> Restart data for interaction containers
   type(container_cache), intent(inout) :: icache

   !> Self-consistent energy
   real(wp), intent(inout) :: energies(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   ...
*/

__device__
void next_scf(
  int &iscf,
  const structure_type &mol,
  const basis_type &bas,
  wavefunction_type wfn,
  sygvd_solver &solver,
  broyden_mixer &mixer,
  const scf_info &info,
  const tb_coulomb &coulomb,
  /*const dispersion_type &dispersion,*/
  /*const container_list &interactions,*/
  const integral_type &ints,
  potential_type &pot,
  coulomb_cache &cache,
  /*container_cache &dcache,*/
  /*container_cache &icache,*/
  float (&energies)[MAX_NAT]
  /*error_type &error*/
);



#endif