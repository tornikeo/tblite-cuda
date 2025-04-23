#ifndef SCF_ITERATOR_H
#define SCF_ITERATOR_H
#include "../basis/type.h"
#include "../limits.h"
#include "../structure.h"
#include "info.h"
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
constexpr int get_mixer_dimension(
  const structure_type &mol, 
  const basis_type &bas,
  const scf_info &info
) {
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

#endif