#ifndef H0_H
#define H0_H
#include "limits.h"
#include "structure.h"
#include "basis/type.h"
#include "integral/type.h"
#include "adjlist.h"
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
   end type tb_hamiltonian*/
typedef struct
{
  /* data */
  float selfenergy[MSHELL][MAX_NELEM];
  float kcn[MSHELL][MAX_NELEM];
  float kq1[MSHELL][MAX_NELEM];
  float kq2[MSHELL][MAX_NELEM];
  float hscale[MSHELL][MSHELL][MAX_NELEM][MAX_NELEM];
  float shpoly[MSHELL][MAX_NELEM];
  float rad[MAX_NELEM];
  float refocc[MSHELL][MAX_NELEM];
  int mshell = MSHELL;
  int mol_nid = MAX_NELEM;
} tb_hamiltonian;

typedef struct
{
  /* data */
  float kshell[2][2];
  float kpair[MAX_NELEM][MAX_NELEM];

} gfn2_h0spec;


/* 
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

__device__
void new_hamiltonian(
  tb_hamiltonian &self,
  const structure_type &mol,
  const basis_type &bas,
  const gfn2_h0spec &spec
);

/*  subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Occupation number
      real(wp), intent(out) :: nocc
      !> Reference occupation for each atom
      real(wp), intent(out) :: n0at(:)
      !> Reference occupation for each shell
      real(wp), intent(out) :: n0sh(:)
*/

__device__
void get_occupation(
  const structure_type &mol,
  const basis_type &bas,
  const tb_hamiltonian &h0,
  float &nocc,
  float (&n0at)[MAX_NAT],
  float (&n0sh)[MAX_NSH]
);


/* 
   subroutine get_selfenergy(h0, id, ish_at, nshell, cn, qat, selfenergy, dsedcn, dsedq)
      type(tb_hamiltonian), intent(in) :: h0
      integer, intent(in) :: id(:)
      integer, intent(in) :: ish_at(:)
      integer, intent(in) :: nshell(:)
      real(wp), intent(in), optional :: cn(:)
      real(wp), intent(in), optional :: qat(:)
      real(wp), intent(out) :: selfenergy(:)
      real(wp), intent(out), optional :: dsedcn(:)
      real(wp), intent(out), optional :: dsedq(:)

      integer :: iat, izp, ish, ii

      selfenergy(:) = 0.0_wp
      if (present(dsedcn)) dsedcn(:) = 0.0_wp
      if (present(dsedq)) dsedq(:) = 0.0_wp
      do iat = 1, size(id)
         izp = id(iat)
         ii = ish_at(iat)
         do ish = 1, nshell(izp)
            selfenergy(ii+ish) = h0%selfenergy(ish, izp)
         end do
      end do
      if (present(cn)) then
         if (present(dsedcn)) then
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
                  dsedcn(ii+ish) = -h0%kcn(ish, izp)
               end do
            end do
          end if
         end if
      end if
  end subroutine get_selfenergy
*/

__device__
void get_selfenergy(
  const tb_hamiltonian &h0, /* TODO: GO FROM HERE! */
  const int (&id)[MAX_NAT],
  const int (&ish_at)[MAX_NAT],
  const int (&nshell)[MAX_NELEM],
  const float (&cn)[MAX_NAT],
  // const float (&qat)[MAX_NAT],
  float selfenergy[MAX_NSH], // static size
  float dsedcn[MAX_NSH]     // static size
  // float dsedq[MAX_NSH],      // static size
  // int id_size,
  // bool cn_present,
  // bool dsedcn_present,
  // bool dsedq_present
);


/* 
subroutine get_hamiltonian(mol, trans, alist, bas, h0, selfenergy, overlap, dpint, qpint, &
   & hamiltonian)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: trans(:, :)
      !> Neighbour list
      type(adjacency_list), intent(in) :: alist
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Diagonal elememts of the Hamiltonian
      real(wp), intent(in) :: selfenergy(:)
      !> Overlap integral matrix
      real(wp), intent(out) :: overlap(:, :)
      !> Dipole moment integral matrix
      real(wp), intent(out) :: dpint(:, :, :)
      !> Quadrupole moment integral matrix
      real(wp), intent(out) :: qpint(:, :, :)
      !> Effective Hamiltonian
      real(wp), intent(out) :: hamiltonian(:, :)

*/

/* 
typedef struct 
{
  float hamiltonian[MAX_NAO][MAX_NAO];
  float overlap[MAX_NAO][MAX_NAO];
  float dipole[MAX_NAO][MAX_NAO][3];  // notice the inverted axes
  float quadrupole[MAX_NAO][MAX_NAO][6]; // Notice the inverted axes
} integral_type;
*/

__device__
void get_hamiltonian(
  const structure_type &mol,
  /* const float (&trans)[MAX_TRANS][3],*/ // Unimplemented
  const adjacency_list &alist,
  const basis_type &bas,
  const tb_hamiltonian &h0,
  const float (&selfenergy)[MAX_NSH],
  float (&overlap)[MAX_NAO][MAX_NAO],
  float (&dpint)[MAX_NAO][MAX_NAO][3],
  float (&qpint)[MAX_NAO][MAX_NAO][6],
  float (&hamiltonian)[MAX_NAO][MAX_NAO]
);


#endif