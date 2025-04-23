#ifndef COULOMB_H
#define COULOMB_H
#include "limits.h"
#include "structure.h"
#include "potential.h"
#include "wavefunction/type.h"

__device__
void get_mrad(
  const structure_type &mol,
  const float shift,
  const float kexp,
  const float rmax,
  const float (&rad)[MAX_NELEM],
  const float (&valence_cn)[MAX_NELEM],
  const float (&cn)[MAX_NAT],
  float (&mrad)[MAX_NAT],
  float (&dmrdcn)[MAX_NAT]
);

/*   type :: coulomb_cache
      real(wp) :: alpha
      type(wignerseitz_cell) :: wsc
      real(wp), allocatable :: amat(:, :)
      real(wp), allocatable :: vvec(:)

      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: dcndr(:, :, :)
      real(wp), allocatable :: dcndL(:, :, :)
      real(wp), allocatable :: mrad(:)
      real(wp), allocatable :: dmrdcn(:)

      real(wp), allocatable :: amat_sd(:, :, :)
      real(wp), allocatable :: amat_dd(:, :, :, :)
      real(wp), allocatable :: amat_sq(:, :, :)
   contains
      procedure :: update
   end type coulomb_cache
*/
typedef struct {
  float alpha;
  //  cell is unused
  // other vars are not used for now
  /*
   if (.not.allocated(ptr%mrad)) then
      allocate(ptr%mrad(mol%nat))
   end if
   if (.not.allocated(ptr%dmrdcn)) then
      allocate(ptr%dmrdcn(mol%nat))
   end if

   if (.not.allocated(ptr%amat_sd)) then
      allocate(ptr%amat_sd(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%amat_dd)) then
      allocate(ptr%amat_dd(3, mol%nat, 3, mol%nat))
   end if
   if (.not.allocated(ptr%amat_sq)) then
      allocate(ptr%amat_sq(6, mol%nat, mol%nat))
   end if

   if (.not.allocated(ptr%cn)) then
      allocate(ptr%cn(mol%nat))
   end if
   if (.not.allocated(ptr%dcndr)) then
      allocate(ptr%dcndr(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%dcndL)) then
      allocate(ptr%dcndL(3, 3, mol%nat))
   end if
  */
  float mrad[MAX_NAT];
  float dmrdcn[MAX_NAT];
  float amat_sd[MAX_NAT][MAX_NAT][3];
  float amat_dd[MAX_NAT][3][MAX_NAT][3];
  float amat_sq[MAX_NAT][MAX_NAT][6];
  float cn[MAX_NAT];
  float dcndr[MAX_NAT][MAX_NAT][3];
  float dcndL[MAX_NAT][3][3];
  float vvec[MSHELL];
  float amat[MAX_NSH][MAX_NSH];
} coulomb_cache;




class coulomb_charge_type 
{
  // public:
  // __device__ void update(const structure_type &mol, container_cache &cache);
};

typedef struct 
{
  float cutoff;
  float rcov[MAX_NELEM];
} gfn_ncoord_type;

/*
type, extends(coulomb_type) :: damped_multipole
      !> Damping function for inverse quadratic contributions
      real(wp) :: kdmp3 = 0.0_wp
      !> Damping function for inverse cubic contributions
      real(wp) :: kdmp5 = 0.0_wp
      !> Kernel for on-site dipole exchange-correlation
      real(wp), allocatable :: dkernel(:)
      !> Kernel for on-site quadrupolar exchange-correlation
      real(wp), allocatable :: qkernel(:)

      !> Shift for the generation of the multipolar damping radii
      real(wp) :: shift = 0.0_wp
      !> Exponent for the generation of the multipolar damping radii
      real(wp) :: kexp = 0.0_wp
      !> Maximum radius for the multipolar damping radii
      real(wp) :: rmax = 0.0_wp
      !> Base radii for the multipolar damping radii
      real(wp), allocatable :: rad(:)
      !> Valence coordination number
      real(wp), allocatable :: valence_cn(:)

      !> Coordination number container for multipolar damping radii
      type(gfn_ncoord_type), allocatable :: ncoord
   contains
   */

class damped_multipole
{ 
  public:
  float kdmp3 = 0;
  float kdmp5 = 0;
  float shift = 0;
  float kexp = 0;
  float rmax = 0;
  float rad[MAX_NELEM];
  float valence_cn[MAX_NELEM];
  float dkernel[MAX_NAT];
  float qkernel[MAX_NAT];

  gfn_ncoord_type ncoord;
  __device__ void update(const structure_type &mol, coulomb_cache &cache) const;
  __device__ void get_multipole_matrix(
    // const damped_multipole &self,
    const structure_type &mol,
    coulomb_cache &cache,
    float (&amat_sd)[MAX_NAT][MAX_NAT][3],
    float (&amat_dd)[MAX_NAT][3][MAX_NAT][3],
    float (&amat_sq)[MAX_NAT][MAX_NAT][6]
  ) const;
};

// class onsite_thirdorder
// {

// };

class effective_coulomb /* TODO: Good candidate for upgrading to a class */
{ 
  public:
  float hubbard[MSHELL][MSHELL][MAX_NELEM][MAX_NELEM];
  int nshell[MAX_NAT];
  int offset[MAX_NSH];
  float gexp;
  float rcut;
  __device__ void update(const structure_type &mol, coulomb_cache &cache) const;
} ;


class tb_coulomb 
{
  public:
    effective_coulomb es2;
    damped_multipole aes2;
    __device__ void update(const structure_type &mol, coulomb_cache &cache) const;
    __device__ void get_potential(
      const structure_type &mol, 
      coulomb_cache &cache,
      const wavefunction_type &wfn,
      potential_type &pot
    ) const;
  /* onsite_thirdorder es3; */ // Unused
};
__device__
void update(
  const tb_coulomb &self, 
  const structure_type &mol,
  coulomb_cache &cache
);

__device__
void ncoord_get_cn(
  const gfn_ncoord_type &self,
  const structure_type &mol,
  float (&cn)[MAX_NAT],
  float (&dcndr)[MAX_NAT][MAX_NAT][3],
  float (&dcndL)[MAX_NAT][3][3]
);

#endif