#include "slater.h"
#include "../limits.h"
#include "type.h"

/*

!> Expand Slater function in primitive gaussian functions
pure subroutine slater_to_gauss_array(ng, n, l, zeta, alpha, coeff, norm, info)
   !> Number of Gaussian functions for the expansion
   integer, intent(in) :: ng
   !> Principal quantum number of shell
   integer, intent(in) :: n
   !> Azimudal quantum number of shell
   integer, intent(in) :: l
   !> Exponent of Slater function to expand
   real(wp), intent(in) :: zeta
   !> Exponents of primitive gaussian functions
   real(wp), intent(out) :: alpha(:)
   !> Contraction coefficients of primitive gaussians, can contain normalization
   real(wp), intent(out) :: coeff(:)
   !> Include normalization in contraction coefficients
   logical, intent(in) :: norm
   !> Status of the expansion, returns zero for success or the position of the
   !> faulty dummy argument
   integer, intent(out) :: info

   ! Storage location of the STO-NG coefficients and exponents
   integer :: ityp

   ! Basic checks first, we cannot support principal QN of six or higher
   ! also you should not violate l ∊ [n-1, n-2, ..., 1, 0]
   if ((n > 5).or.(n <= l)) then
      if (.not.(n == 6 .and. ng == 6)) then
         info = 2
         return
      end if
   endif

   ! Don't allow negative exponents in the first place
   if (zeta <= 0.0_wp) then
      info = 4
      return
   endif

   ! we have to use a little hack here,
   ! if you pass n and l correctly, everything is fine
   ! ityp: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
   !    n: 1 2 3 4 5 2 3 4 5  3  4  5  4  5  5
   !    l: 0 0 0 0 0 1 1 1 1  2  2  2  3  3  4
   select case(l) ! integer hack:
   case(0); ityp = n    ! s
   case(1); ityp = 4+n  ! p
   case(2); ityp = 7+n  ! d
   case(3); ityp = 9+n  ! f
   case(4); ityp = 10+n ! g
   case default ! I'm sorry, no h-functions for you today
      info = 3
      return
   end select

   select case(ng)
   case default ! currently we cannot go beyond 6 primitives
      info = 1
      return
   case(1)
      alpha(1) = pAlpha1(ityp) * zeta**2
      coeff(1) = 1.0_wp
   case(2)
      alpha(:ng) = pAlpha2(:, ityp) * zeta**2
      coeff(:ng) = pCoeff2(:, ityp)
   case(3)
      alpha(:ng) = pAlpha3(:, ityp) * zeta**2
      coeff(:ng) = pCoeff3(:, ityp)
   case(4)
      alpha(:ng) = pAlpha4(:, ityp) * zeta**2
      coeff(:ng) = pCoeff4(:, ityp)
   case(5)
      alpha(:ng) = pAlpha5(:, ityp) * zeta**2
      coeff(:ng) = pCoeff5(:, ityp)
   case(6)
      if (n == 6) then
         if (l == 0) then
            alpha(:ng) = pAlpha6s(:) * zeta**2
            coeff(:ng) = pCoeff6s(:)
         else if (l == 1) then
            alpha(:ng) = pAlpha6p(:) * zeta**2
            coeff(:ng) = pCoeff6p(:)
         else
            info = 2
            return
         end if
      else
         alpha(:ng) = pAlpha6(:, ityp) * zeta**2
         coeff(:ng) = pCoeff6(:, ityp)
      end if
   end select

   ! normalize the gaussian if requested
   ! <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
   ! N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
   ! N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   if (norm) then
      coeff(:ng) = coeff(:ng) * (top*alpha(:ng))**0.75_wp &
         & * sqrt(4*alpha(:ng))**l / sqrt(dfactorial(l+1))
   endif

   ! success
   info = 0

end subroutine slater_to_gauss_array
*/

// void slater_to_gauss(
//   const int ng, const int n, const int l, const float zeta, 
//   float (&alpha)[MAXG], 
//   float (&coeff)[MAXG],
//   const bool norm,
//   int &info
// )
// {
//   // Basic checks
//   if (n > 5 || n <= l) {
//     if (!(n == 6 && ng == 6)) {
//       info = 2;
//       return;
//     }
//   }

//   if (zeta <= 0.0f) {
//     info = 4;
//     return;
//   }

//   // Determine ityp based on n and l
//   int ityp;
//   switch (l) {
//     case 0: ityp = n; break;       // s
//     case 1: ityp = 4 + n; break;   // p
//     case 2: ityp = 7 + n; break;   // d
//     case 3: ityp = 9 + n; break;   // f
//     case 4: ityp = 10 + n; break;  // g
//     default:
//       info = 3;
//       return;
//   }

//   // Handle ng cases
//   switch (ng) {
//     default:
//       info = 1;
//       return;
//     case 1:
//       alpha[0] = pAlpha1[ityp] * zeta * zeta;
//       coeff[0] = 1.0f;
//       break;
//     case 2:
//       for (int i = 0; i < ng; ++i) {
//         alpha[i] = pAlpha2[i][ityp] * zeta * zeta;
//         coeff[i] = pCoeff2[i][ityp];
//       }
//       break;
//     case 3:
//       for (int i = 0; i < ng; ++i) {
//         alpha[i] = pAlpha3[i][ityp] * zeta * zeta;
//         coeff[i] = pCoeff3[i][ityp];
//       }
//       break;
//     case 4:
//       for (int i = 0; i < ng; ++i) {
//         alpha[i] = pAlpha4[i][ityp] * zeta * zeta;
//         coeff[i] = pCoeff4[i][ityp];
//       }
//       break;
//     case 5:
//       for (int i = 0; i < ng; ++i) {
//         alpha[i] = pAlpha5[i][ityp] * zeta * zeta;
//         coeff[i] = pCoeff5[i][ityp];
//       }
//       break;
//     case 6:
//       if (n == 6) {
//         if (l == 0) {
//           for (int i = 0; i < ng; ++i) {
//             alpha[i] = pAlpha6s[i] * zeta * zeta;
//             coeff[i] = pCoeff6s[i];
//           }
//         } else if (l == 1) {
//           for (int i = 0; i < ng; ++i) {
//             alpha[i] = pAlpha6p[i] * zeta * zeta;
//             coeff[i] = pCoeff6p[i];
//           }
//         } else {
//           info = 2;
//           return;
//         }
//       } else {
//         for (int i = 0; i < ng; ++i) {
//           alpha[i] = pAlpha6[i][ityp] * zeta * zeta;
//           coeff[i] = pCoeff6[i][ityp];
//         }
//       }
//       break;
//   }

//   // Normalize if requested
//   if (norm) {
//     for (int i = 0; i < ng; ++i) {
//       coeff[i] *= powf(top * alpha[i], 0.75f) *
//                   powf(sqrtf(4 * alpha[i]), l) /
//                   sqrtf(dfactorial[l + 1]);
//     }
//   }

//   // Success
//   info = 0;
// }