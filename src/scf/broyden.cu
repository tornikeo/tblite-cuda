#include "broyden.h"

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

/*
subroutine set_1d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), intent(in) :: qvec(:)
   self%q_in(self%iset+1:self%iset+size(qvec)) = qvec
   self%iset = self%iset + size(qvec)
end subroutine set_1d*/

__device__
void broyden_mixer::set(const float *qvec, int size)
{
  for (int i = 0; i < size; ++i) {
    q_in[iset + i] = qvec[i];
  }
  iset += size;
}
/*

subroutine diff_1d(self, qvec)
   class(broyden_mixer), intent(inout) :: self
   real(wp), intent(in) :: qvec(:)
   self%dq(self%idif+1:self%idif+size(qvec)) = qvec &
      & - self%q_in(self%idif+1:self%idif+size(qvec))
   self%idif = self%idif + size(qvec)
end subroutine diff_1d
*/
__device__
void broyden_mixer::diff(
  const float *qsh, int size
)
{
  for (int i = idif; i < idif + size; i++)
  {
    dq[i] = qsh[i] - q_in[i];
  }
  idif += size;
}
