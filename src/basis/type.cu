#include "type.h"
#include "limits.h"

/*

pure function clip(val, min_val, max_val) result(res)
   real(wp), intent(in) :: val, min_val, max_val
   real(wp) :: res
   res = min(max(val, min_val), max_val)
end function clip
*/
__device__
float clip(const float val, const float min_val, const float max_val)
{
  return fminf(fmaxf(val, min_val), max_val);
}
/*
pure function integral_cutoff(acc) result(intcut)
   !> Accuracy for the integral cutoff
   real(wp), intent(in) :: acc
   !> Integral cutoff
   real(wp) :: intcut

   real(wp), parameter :: min_intcut = 5.0_wp, max_intcut = 25.0_wp, &
      & max_acc = 1.0e-4_wp, min_acc = 1.0e+3_wp

   intcut = clip(max_intcut - 10*log10(clip(acc, min_acc, max_acc)), min_intcut, max_intcut)
end function integral_cutoff

*/

__device__
float integral_cutoff(const float acc)
{
  const float min_intcut = 5.0f;
  const float max_intcut = 25.0f;
  const float max_acc = 1.0e-4f;
  const float min_acc = 1.0e+3f;

  return clip(max_intcut - 10.0f * log10f(clip(acc, min_acc, max_acc)), min_intcut, max_intcut);
}

/* pure function get_cutoff(self, acc) result(cutoff)
   !> Instance of the basis set data
   type(basis_type), intent(in) :: self
   !> Accuracy for the integral cutoff
   real(wp), intent(in), optional :: acc
   !> Required realspace cutoff
   real(wp) :: cutoff

   real(wp) :: intcut
   real(wp), parameter :: max_cutoff = 40.0_wp

   if (present(acc)) then
      intcut = integral_cutoff(acc)
   else
      intcut = self%intcut
   end if
   ! ai * aj * cutoff2 / (ai + aj) == intcut
   cutoff = min(sqrt(2.0_wp*intcut/self%min_alpha), max_cutoff)

end function get_cutoff

*/

__device__
float get_cutoff(const basis_type &self, const float acc)
{
  const float max_cutoff = 40.0f;
  float intcut;

  if (acc >= 0.0f) {
    intcut = integral_cutoff(acc);
  } else {
    intcut = self.intcut;
  }

  return fminf(sqrtf(2.0f * intcut / self.min_alpha), max_cutoff);
}