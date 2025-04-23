#include <stdio.h>
#include "iterators.h"
#include <cassert>

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
  coulomb_cache &cache, /* NOTE: Not an abstract container_cache since we don't mess with the polymorphism for now */
  /*container_cache &dcache,*/
  /*container_cache &icache,*/
  float (&energies)[MAX_NAT]
  /*error_type &error*/
) {
  printf("NEXT SCF %i\n", iscf);
  float eao[MAX_NAO] = {0};
  float ts = 0;
  if (iscf > 0)
  {
    // next(mixer);
    // get_mixer(mixer, bas, wfn, info); // Update wfn from mixer, pretty much
    // assert(false && "Unimplemented");
  }
  pot.reset(); /* Dipping toes in class member functions */
  iscf++;

  /* if (present(coulomb)) then */
    // call coulomb%get_potential(mol, cache, wfn, pot)
    coulomb.get_potential(mol, cache, wfn, pot);
  /*end if*/
  /*if (present(dispersion)) then
    call dispersion%get_potential(mol, dcache, wfn, pot)
  end if*/
  /*if (present(interactions)) then
    call interactions%get_potential(mol, icache, wfn, pot)
  end if*/
}