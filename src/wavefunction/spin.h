#ifndef WAVEFUNCTION_SPIN_H
#define WAVEFUNCTION_SPIN_H
#include "../limits.h"
/*
subroutine updown_to_magnet_2(x)

   !> Data in up-down representation on entry
   !> transformed to charge-magnetization representation on exit
   real(wp), intent(inout) :: x(:, :)

   if (size(x, 2) == 2) then
      x(:, 1) = x(:, 1) + x(:, 2)
      x(:, 2) = x(:, 1) - 2.0_wp * x(:, 2)
   end if
end subroutine updown_to_magnet_2
*/

__device__
void updown_to_magnet(
  float (&x)[MAX_NSPIN][MAX_NSH]
);

__device__
void updown_to_magnet(
  float (&x)[MAX_NSPIN][MAX_NAT][3]
);

#endif