#ifndef INTEGRAL_MUTLIPOLE_H
#define INTEGRAL_MUTLIPOLE_H
#include "../limits.h"
#include "../basis/type.h"
#include "overlap.h"
/*

elemental function overlap_1d(moment, alpha) result(overlap)
   integer, intent(in) :: moment
   real(wp), intent(in) :: alpha
   real(wp) :: overlap
   real(wp), parameter :: dfactorial(0:7) = & ! see OEIS A001147
      & [1._wp,1._wp,3._wp,15._wp,105._wp,945._wp,10395._wp,135135._wp]

   if (modulo(moment, 2) == 0) then
      overlap = (0.5_wp/alpha)**(moment/2) * dfactorial(moment/2)
   else
      overlap = 0.0_wp
   end if
end function overlap_1d

*/
// ! x (+1), y (-1), z (0) in [-1, 0, 1] sorting
// integer, parameter :: lx(3, 84) = reshape([&
//   & 0, &
//   & 0,0,1, &
//   & 2,0,0,1,1,0, &
//   & 3,0,0,2,2,1,0,1,0,1, &
//   & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
//   & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
//   & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2, &
//   & 0, &
//   & 1,0,0, &
//   & 0,2,0,1,0,1, &
//   & 0,3,0,1,0,2,2,0,1,1, &
//   & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
//   & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
//   & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2, &
//   & 0, &
//   & 0,1,0, &
//   & 0,0,2,0,1,1, &
//   & 0,0,3,0,1,0,1,2,2,1, &
//   & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
//   & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
//   & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2], &
//   & shape(lx), order=[2, 1])

/*
pure subroutine multipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6), pre, tr
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3d(6, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3d(:, :, :) = 0.0_wp
   q3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 2
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, quad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3d(:, mlj, mli) = d3d(:, mlj, mli) + cc*dip
               q3d(:, mlj, mli) = q3d(:, mlj, mli) + cc*quad
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap)
   call transform1(cgtoj%ang, cgtoi%ang, d3d, dpint)
   call transform1(cgtoj%ang, cgtoi%ang, q3d, qpint)

   ! remove trace from quadrupole integrals (transfrom to spherical harmonics and back)
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do

end subroutine multipole_cgto

*/


/*
pure subroutine multipole_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint, &
      & doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpinti(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpinti(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpintj(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpintj(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6)
   real(wp) :: pre, grad(3), ddip(3, 3), dquad(3, 6), tr, dtr(3)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3dj(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3dj(6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3di(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3di(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3dj(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3dj(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   
   s3d(:, :) = 0.0_wp
   d3dj(:, :, :) = 0.0_wp
   q3dj(:, :, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp
   dd3di(:, :, :, :) = 0.0_wp
   dq3di(:, :, :, :) = 0.0_wp
   dd3dj(:, :, :, :) = 0.0_wp
   dq3dj(:, :, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 3
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_grad_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, quad, grad, ddip, dquad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3dj(:, mlj, mli) = d3dj(:, mlj, mli) + cc*dip
               q3dj(:, mlj, mli) = q3dj(:, mlj, mli) + cc*quad
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
               dd3dj(:, :, mlj, mli) = dd3dj(:, :, mlj, mli) + cc*ddip
               dq3dj(:, :, mlj, mli) = dq3dj(:, :, mlj, mli) + cc*dquad
            end do
         end do
      end do
   end do

   do mli = 1, mlao(cgtoi%ang)
      do mlj = 1, mlao(cgtoj%ang)
         call shift_operator(vec, s3d(mlj, mli), d3dj(:, mlj, mli), q3dj(:, mlj, mli), &
            & ds3d(:, mlj, mli), dd3dj(:, :, mlj, mli), dq3dj(:, :, mlj, mli), &
            & dd3di(:, :, mlj, mli), dq3di(:, :, mlj, mli))
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap)
   call transform1(cgtoj%ang, cgtoi%ang, d3dj, dpint)
   call transform1(cgtoj%ang, cgtoi%ang, q3dj, qpint)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap)
   call transform2(cgtoj%ang, cgtoi%ang, dd3dj, ddpintj)
   call transform2(cgtoj%ang, cgtoi%ang, dq3dj, dqpintj)
   call transform2(cgtoj%ang, cgtoi%ang, dd3di, ddpinti)
   call transform2(cgtoj%ang, cgtoi%ang, dq3di, dqpinti)

   ! remove trace from quadrupole integrals (transfrom to spherical harmonics and back)
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do

   ! remove trace from quadrupole integral derivatives
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         dtr = dqpinti(:, 1, mlj, mli) + dqpinti(:, 3, mlj, mli) + dqpinti(:, 6, mlj, mli)
         dqpinti(:, 1, mlj, mli) = 1.5_wp * dqpinti(:, 1, mlj, mli) - 0.5_wp * dtr
         dqpinti(:, 2, mlj, mli) = 1.5_wp * dqpinti(:, 2, mlj, mli)
         dqpinti(:, 3, mlj, mli) = 1.5_wp * dqpinti(:, 3, mlj, mli) - 0.5_wp * dtr
         dqpinti(:, 4, mlj, mli) = 1.5_wp * dqpinti(:, 4, mlj, mli)
         dqpinti(:, 5, mlj, mli) = 1.5_wp * dqpinti(:, 5, mlj, mli)
         dqpinti(:, 6, mlj, mli) = 1.5_wp * dqpinti(:, 6, mlj, mli) - 0.5_wp * dtr
      end do
   end do

   ! remove trace from quadrupole integral derivatives
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         dtr = dqpintj(:, 1, mlj, mli) + dqpintj(:, 3, mlj, mli) + dqpintj(:, 6, mlj, mli)
         dqpintj(:, 1, mlj, mli) = 1.5_wp * dqpintj(:, 1, mlj, mli) - 0.5_wp * dtr
         dqpintj(:, 2, mlj, mli) = 1.5_wp * dqpintj(:, 2, mlj, mli)
         dqpintj(:, 3, mlj, mli) = 1.5_wp * dqpintj(:, 3, mlj, mli) - 0.5_wp * dtr
         dqpintj(:, 4, mlj, mli) = 1.5_wp * dqpintj(:, 4, mlj, mli)
         dqpintj(:, 5, mlj, mli) = 1.5_wp * dqpintj(:, 5, mlj, mli)
         dqpintj(:, 6, mlj, mli) = 1.5_wp * dqpintj(:, 6, mlj, mli) - 0.5_wp * dtr
      end do
   end do
end subroutine multipole_grad_cgto
*/

__device__
void multipole_cgto(
  const cgto_type &cgtoi,
  const cgto_type &cgtoj,
  const float r2,
  const float (&vec)[3],
  const float intcut,
  float (&overlap)[msao(MAXL)][msao(MAXL)],
  float (&dpint)[msao(MAXL)][msao(MAXL)][3],
  float (&qpint)[msao(MAXL)][msao(MAXL)][6]
);

#endif 