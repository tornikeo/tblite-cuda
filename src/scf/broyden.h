#ifndef SCF_BROYDEN_H
#define SCF_BROYDEN_H
#include "../limits.h"
/*
   type :: broyden_mixer
      integer :: ndim
      integer :: memory
      integer :: iter
      integer :: iset
      integer :: idif
      integer :: iget
      real(wp) :: damp
      real(wp), allocatable :: df(:, :)
      real(wp), allocatable :: u(:, :)
      real(wp), allocatable :: a(:, :)
      real(wp), allocatable :: dq(:)
      real(wp), allocatable :: dqlast(:)
      real(wp), allocatable :: qlast_in(:)
      real(wp), allocatable :: omega(:)
      real(wp), allocatable :: q_in(:)
   contains
*/


/*
subroutine new_broyden(self, memory, ndim, damp)
   type(broyden_mixer), intent(out) :: self
   integer, intent(in) :: memory
   integer, intent(in) :: ndim
   real(wp), intent(in) :: damp

   self%ndim = ndim
   self%memory = memory
   self%iter = 0
   self%iset = 0
   self%idif = 0
   self%iget = 0
   self%damp = damp
   allocate(self%df(ndim, memory))
   allocate(self%u(ndim, memory))
   allocate(self%a(memory, memory))
   allocate(self%dq(ndim))
   allocate(self%dqlast(ndim))
   allocate(self%qlast_in(ndim))
   allocate(self%omega(memory))
   allocate(self%q_in(ndim))
end subroutine new_broyden
*/

typedef struct{
  int ndim = BROYDEN_NDIM;
  int memory = MAX_ITER_DEFAULT;
  int iter = 0;
  int iset = 0;
  int idif = 0;
  int iget = 0;
  float damp = 0.0f;
  float df[BROYDEN_NDIM][MAX_ITER_DEFAULT] = {0}; /* TODO: Prio med. These are some of the 
  largest local memory allocations here. We need to make these smaller */
  float u[BROYDEN_NDIM][MAX_ITER_DEFAULT] = {0};
  float a[MAX_ITER_DEFAULT][MAX_ITER_DEFAULT] = {0};
  float dq[BROYDEN_NDIM] = {0};
  float dqlast[BROYDEN_NDIM] = {0};
  float qlast_in[BROYDEN_NDIM] = {0};
  float omega[MAX_ITER_DEFAULT] = {0};
  float q_in[BROYDEN_NDIM] = {0};
} broyden_mixer;


__device__
void new_broyden(
  broyden_mixer &self,
  const int memory,
  const int ndim,
  const float damp
) {
  self.ndim = ndim;
  self.memory = memory;
  self.iter = 0;
  self.iset = 0;
  self.idif = 0;
  self.iget = 0;
  self.damp = damp;

  // Allocate memory for arrays
  // self.df = new float[ndim * memory];
  // self.u = new float[ndim * memory];
  // self.a = new float[memory * memory];
  // self.dq = new float[ndim];
  // self.dqlast = new float[ndim];
  // self.qlast_in = new float[ndim];
  // self.omega = new float[memory];
  // self.q_in = new float[ndim];

  // // Initialize allocated memory to zero
  // for (int i = 0; i < ndim * memory; ++i) {
  //   self.df[i] = 0.0f;
  //   self.u[i] = 0.0f;
  // }
  // for (int i = 0; i < memory * memory; ++i) {
  //   self.a[i] = 0.0f;
  // }
  // for (int i = 0; i < ndim; ++i) {
  //   self.dq[i] = 0.0f;
  //   self.dqlast[i] = 0.0f;
  //   self.qlast_in[i] = 0.0f;
  //   self.q_in[i] = 0.0f;
  // }
  // for (int i = 0; i < memory; ++i) {
  //   self.omega[i] = 0.0f;
  // }
}

#endif