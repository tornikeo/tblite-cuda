#include "xtb.h"
#include "gfn2.h"
#include "limits.h"
#include "utils/array.h"
#include "structure.h"
#include "basis/slater.h"
#include <stdio.h>

void init(structure_type &mol)
{

  int nat =            9 ;
  int nid =            4 ;
  int nbd =            0 ;
  int id[] = {
 1, 1, 2, 3, 4, 4, 4, 4, 4 };
  int num[] = {
 6, 17, 8, 1 };
//   const char* sym[] = {
//  "C", "Cl", "O", "H" };
  float xyz[9][3] = {
   {0.98819706, -0.79388595, 1.8932948 },
   {-1.1710542, 0.74677960E-1, 3.3566231 },
   {1.6169206, -3.8529364, 3.2129892 },
   {-1.2120192, 2.6964878, 3.6464740 },
   {-1.6371701, 3.4899489, 2.1451411 },
   {-2.9580292, -0.57318774, 2.5348828 },
   {-0.97500824, -0.79775967, 5.2966186 },
   {0.40421389, -1.4110043, -0.10239631 },
   {2.7573567, 0.39135773, 1.9942774 },
  };
  int uhf =            0 ;
  double charge =    1.0000000000000000;

  mol.nat = nat;
  mol.nid = nid;
  mol.nbd = nbd;

  set(mol.id, id);
  set(mol.num, num);
  // set(mol.sym, sym);
  set(mol.xyz, xyz);
  mol.uhf = uhf;
  mol.charge = charge;
  mol.periodic = false;
}

void init(wavefunction_type &wfn)
{
  
}

/*subroutine add_basis(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: isp, izp, ish, stat, ng
   integer, allocatable :: nsh_id(:)
   type(cgto_type), allocatable :: cgto(:, :)

   nsh_id = nshell(mol%num)
   allocate(cgto(maxval(nsh_id), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nsh_id(isp)
         ng = number_of_primitives(ish, izp)
         call slater_to_gauss(ng, principal_quantum_number(ish, izp), ang_shell(ish, izp), &
            & slater_exponent(ish, izp), cgto(ish, isp), .true., stat)
      end do
   end do

   call new_basis(calc%bas, mol, nsh_id, cgto, 1.0_wp)

end subroutine add_basis
*/



void add_basis(xtb_calculator &calc, structure_type &mol)
{
  int nsh_id[MAX_NELEM];
  cgto_type cgto[MAX_NSH][MAX_NELEM];

  for (int isp = 0; isp < mol.nid; isp++)
  {
    int izp = mol.num[isp];
    nsh_id[isp] = nshell[izp];

    for (int ish = 0; ish < nsh_id[isp]; ish++)
    {
      int ng = number_of_primitives[ish][izp];
      cgto[ish][isp];
      int stat = 0;
      slater_to_gauss(ng, principal_quantum_number[ish][izp], ang_shell[ish][izp],
                      slater_exponent[ish][izp], 
                      cgto[ish][isp].alpha, cgto[ish][isp].coeff, 
                      true, stat);
      assert(stat != 0); // add basis failed
    }
  }

  // new_basis(calc.bas, mol, nsh_id, cgto, 1.0);
}

void init(xtb_calculator &calc, structure_type &mol)
{
  add_basis(calc, mol);
  // add_ncoord(calc, mol);
  // add_hamiltonian(calc, mol);
  // add_repulsion(calc, mol);
  // add_dispersion(calc, mol);
  // add_coulomb(calc, mol);
}

int main()
{
  structure_type mol;
  init(mol);

  xtb_calculator calc;
  init(calc, mol);
//   int ang =            0 ;
//  int nprim =            4 ;
//  double alpha[] = {
//   51.049363, 8.7911227, 0.70640422, 0.26922813, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000 };
//  double coeff[] = {
//   -0.16312063, -0.19911025, 0.31882048, 0.12706504, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000 };
  printf("Inside Testing xtb\n");
  
}