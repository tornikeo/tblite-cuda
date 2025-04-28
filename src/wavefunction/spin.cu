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

__device__
void updown_to_magnet(
  float (&x)[MAX_NSPIN][MAX_NSH]
)
{
  if (MAX_NSPIN == 2) { // since fortran indices are reverse of c-indices
    x[0][0] = x[0][0] + x[1][0];
    x[1][0] = x[0][0] - 2.0f * x[1][0];
  }
}


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

__device__
void updown_to_magnet(
  float (&x)[MAX_NSPIN][MAX_NAT][3]
)
{
  if (MAX_NSPIN == 2) { // since fortran indices are reverse of c-indices
    for (int i = 0; i < MAX_NAT; ++i) {
      for (int j = 0; j < 3; ++j) {
        x[0][i][j] = x[0][i][j] + x[1][i][j];
        x[1][i][j] = x[0][i][j] - 2.0f * x[1][i][j];
      }
    }
  }
}