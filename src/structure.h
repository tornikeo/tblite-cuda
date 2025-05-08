#ifndef STRUCTURE_H
#define STRUCTURE_H
#include "limits.h"


constexpr __device__ char* to_symbol(int atomic_num)
{
  switch (atomic_num) {
    case 1:  return "H";
    case 2:  return "He";
    case 3:  return "Li";
    case 4:  return "Be";
    case 5:  return "B";
    case 6:  return "C";
    case 7:  return "N";
    case 8:  return "O";
    case 9:  return "F";
    case 10: return "Ne";
    case 11: return "Na";
    case 12: return "Mg";
    case 13: return "Al";
    case 14: return "Si";
    case 15: return "P";
    case 16: return "S";
    default: return "?";
  }
}

class structure_type
{
  public:
  int nat;
  int nid;
  int nbd;
  bool periodic = false; 
  int id[MAX_NAT];
  int num[MAX_NELEM];
  float xyz[MAX_NAT][3];
  float charge = 0;
  int uhf = 0;
  char* sym[MAX_NAT];

  /*
  subroutine new_structure_num(self, num, xyz, charge, uhf, lattice, periodic, &
      & info, bond)

   !> Instance of the structure representation
   type(structure_type), intent(out) :: self

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Total charge
   real(wp), intent(in), optional :: charge

   !> Number of unpaired electrons
   integer, intent(in), optional :: uhf

   !> Lattice parameters
   real(wp), intent(in), optional :: lattice(:, :)

   !> Periodic directions
   logical, intent(in), optional :: periodic(:)

   !> Vendor specific structure information
   type(structure_info), intent(in), optional :: info

   !> Bond topology of the system
   integer, intent(in), optional :: bond(:, :)

   integer :: ndim, iat
   character(len=symbol_length), allocatable :: sym(:)
*/
  __device__
  inline void new_(
    int num[MAX_NAT],
    float xyz[MAX_NAT][3],
    float charge,
    int uhf,
    float lattice[3][3],
    bool periodic[3],
    int bond[MAX_NAT][MAX_NAT]
  ) {
    int ndim = MAX_NAT;
    for(int iat = 0; iat< ndim; iat++)
    {
      sym[iat] = to_symbol(num[iat]);
    }
  }
};

#endif