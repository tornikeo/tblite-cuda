#include "spin.h"
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



/*

!> Convert an up-down representation into a charge-magnetization representation
subroutine updown_to_magnet_3(x)

   !> Data in up-down representation on entry
   !> transformed to charge-magnetization representation on exit
   real(wp), intent(inout) :: x(:, :, :)

   if (size(x, 3) == 2) then
      x(:, :, 1) = x(:, :, 1) + x(:, :, 2)
      x(:, :, 2) = x(:, :, 1) - 2.0_wp * x(:, :, 2)
   end if
end subroutine updown_to_magnet_3
*/