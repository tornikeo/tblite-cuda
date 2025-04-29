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

template <int Dim>
__device__
void updown_to_magnet(
  float (&x)[MAX_NSPIN][MAX_NAT][Dim]
)
{
  if (MAX_NSPIN == 2) { // since fortran indices are reverse of c-indices
    for (int i = 0; i < MAX_NAT; ++i) {
      for (int j = 0; j < Dim; ++j) {
        x[0][i][j] = x[0][i][j] + x[1][i][j];
        x[1][i][j] = x[0][i][j] - 2.0f * x[1][i][j];
      }
    }
  }
}

#endif