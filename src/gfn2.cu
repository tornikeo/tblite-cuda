#include "gfn2.h"
#include "basis/slater.h"
#include "./limits.h"
#include <cmath>
#include <cstdlib>

/*!> Expand Slater function in primitive gaussian functions
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
#include "gfn2.h"
#include <cmath>
#include <cstdlib>


end subroutine slater_to_gauss_array*/

// Dummy declarations for global arrays and functions. These must be defined elsewhere.
// extern double pAlpha1[];                    // Indexed by ityp
// extern double** pAlpha2;                    // [ng][ityp]
// extern double** pAlpha3;
// extern double** pAlpha4;
// extern double** pAlpha5;
// extern double** pAlpha6;                    // Used when n != 6, or for special cases below
// extern double* pAlpha6s;                    // For n == 6, l == 0
// extern double* pAlpha6p;                    // For n == 6, l == 1

// extern double** pCoeff2;
// extern double** pCoeff3;
// extern double** pCoeff4;
// extern double** pCoeff5;
// extern double** pCoeff6;
// extern double* pCoeff6s;
// extern double* pCoeff6p;
// extern double top;                          // Some constant needed for normalization

// Dummy factorial function: computes the double factorial-like term dfactorial(n)
// double dfactorial(int n) {
//   double result = 1.0;
//   for (int i = 1; i <= n; i++) {
//     result *= (2 * i - 1);
//   }
//   return result;
// }

void slater_to_gauss_array(
  int ng, 
  int n, 
  int l, 
  double zeta, 
  double (&alpha)[MAXG], // output array of length >= ng
  double (&coeff)[MAXG], // output array of length >= ng
  bool norm, 
  int &info)
{
  // Basic input checks
  if ((n > 5 || n <= l) && !(n == 6 && ng == 6)) {
    info = 2;
    return;
  }
  if (zeta <= 0.0) {
    info = 4;
    return;
  }

  // Determine 'ityp' based on l
  int ityp = 0;
  switch(l) {
    case 0: ityp = n;      break; // s
    case 1: ityp = 4 + n;  break; // p
    case 2: ityp = 7 + n;  break; // d
    case 3: ityp = 9 + n;  break; // f
    case 4: ityp = 10 + n; break; // g
    default:
      info = 3;
      return;
  }

  // Expand Slater function into the Gaussian primitives
  switch(ng) {
    default:
      info = 1;
      return;
    case 1:
      alpha[0] = pAlpha1[ityp] * (zeta * zeta);
      coeff[0] = 1.0;
      break;
    case 2:
      for (int i = 0; i < 2; i++) {
        alpha[i] = pAlpha2[i][ityp] * (zeta * zeta);
        coeff[i] = pCoeff2[i][ityp];
      }
      break;
    case 3:
      for (int i = 0; i < 3; i++) {
        alpha[i] = pAlpha3[i][ityp] * (zeta * zeta);
        coeff[i] = pCoeff3[i][ityp];
      }
      break;
    case 4:
      for (int i = 0; i < 4; i++) {
        alpha[i] = pAlpha4[i][ityp] * (zeta * zeta);
        coeff[i] = pCoeff4[i][ityp];
      }
      break;
    case 5:
      for (int i = 0; i < 5; i++) {
        alpha[i] = pAlpha5[i][ityp] * (zeta * zeta);
        coeff[i] = pCoeff5[i][ityp];
      }
      break;
    case 6:
      if (n == 6) {
        if (l == 0) {
          for (int i = 0; i < 6; i++) {
            alpha[i] = pAlpha6s[i] * (zeta * zeta);
            coeff[i] = pCoeff6s[i];
          }
        }
        else if (l == 1) {
          for (int i = 0; i < 6; i++) {
            alpha[i] = pAlpha6p[i] * (zeta * zeta);
            coeff[i] = pCoeff6p[i];
          }
        }
        else {
          info = 2;
          return;
        }
      }
      else {
        for (int i = 0; i < 6; i++) {
          alpha[i] = pAlpha6[i][ityp] * (zeta * zeta);
          coeff[i] = pCoeff6[i][ityp];
        }
      }
      break;
  }

  // Normalize the Gaussian if requested
  if (norm) {
    for (int i = 0; i < ng; i++) {
      // The normalization factor as given in the original code:
      // coeff *= (top * alpha)^0.75 * sqrt((4*alpha)^l) / sqrt(dfactorial(l+1))
      coeff[i] *= std::pow(top * alpha[i], 0.75) *
            std::sqrt(std::pow(4 * alpha[i], l)) /
            std::sqrt(dfactorial[l + 1]);
    }
  }

  // Success
  info = 0;
}
