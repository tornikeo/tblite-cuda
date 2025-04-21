#include "limits.h"
#include "structure.h"

class coulomb_charge_type 
{
  // public:
  // __device__ void update(const structure_type &mol, container_cache &cache);
};

class damped_multipole
{

};

class onsite_thirdorder
{

};

class tb_coulomb {
  public:
  coulomb_charge_type es2;
  damped_multipole aes2;
  onsite_thirdorder es3;
  __device__ void update(const structure_type &mol);
};
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
  float amat[MAX_NSH][MAX_NSH];

} coulomb_cache;


typedef struct 
{ 
  float hubbard[MSHELL][MSHELL][MAX_NELEM][MAX_NELEM];
  int nshell[MAX_NAT];
  int offset[MAX_NSH];
  float gexp;
  float rcut;
} effective_coulomb;


__device__
void update(const tb_coulomb &self, 
  const structure_type &mol,
  coulomb_cache &cache);