#include <stdio.h>
#include "iterators.h"

__device__
// constexpr 
int get_mixer_dimension(
  const structure_type &mol, 
  const basis_type &bas,
  const scf_info &info
)
{
  int ndim = 0;

  switch (info.charge) {
    case atom_resolved:
      ndim += mol.nat;
      break;
    case shell_resolved:
      ndim += bas.nsh;
      break;
  }

  switch (info.dipole) {
    case atom_resolved:
      ndim += 3 * mol.nat;
      break;
  }

  switch (info.quadrupole) {
    case atom_resolved:
      ndim += 6 * mol.nat;
      break;
  }

  return ndim;
}


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
  tb_coulomb &cache,
  /*container_cache &dcache,*/
  /*container_cache &icache,*/
  float (&energies)[MAX_NAT]
  /*error_type &error*/
) {
  printf("NEXT SCF %i\n", iscf);
}