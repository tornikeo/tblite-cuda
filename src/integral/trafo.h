#ifndef INTEGRAL_TRAFO_H
#define INTEGRAL_TRAFO_H
#include "overlap.h"
#include "../limits.h"

/*
pure subroutine transform0(lj, li, cart, sphr)
   integer, intent(in) :: li
   integer, intent(in) :: lj
   real(wp), intent(in) :: cart(:, :)
   real(wp), intent(out) :: sphr(:, :)

   select case(li)
   case(0, 1)
      select case(lj)
      case(0, 1)
         sphr = cart
      case(2)
         ! sphr = matmul(dtrafo, cart)
         ! sphr(1,1) = cart(3, 1) - .5 * (cart(1,1) + cart(2,:))
         ! sphr(1*col + 1) = cart(3 * col + 1) - .5 * (cart(1,1) + cart(2,:))
         sphr(1, :) = cart(3, :) - 0.5_wp * (cart(1, :) + cart(2, :))
         sphr(2, :) = s3 * cart(5, :)
         sphr(3, :) = s3 * cart(6, :)
         sphr(4, :) = s3_4 * (cart(1, :) - cart(2, :))
         sphr(5, :) = s3 * cart(4, :)
      case(3)
         sphr = matmul(ftrafo, cart)
      case(4)
         sphr = matmul(gtrafo, cart)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(2)
      select case(lj)
      case(0, 1)
         ! sphr = matmul(cart, transpose(dtrafo))
         sphr(:, 1) = cart(:, 3) - 0.5_wp * (cart(:, 1) + cart(:, 2))
         sphr(:, 2) = s3 * cart(:, 5)
         sphr(:, 3) = s3 * cart(:, 6)
         sphr(:, 4) = s3_4 * (cart(:, 1) - cart(:, 2))
         sphr(:, 5) = s3 * cart(:, 4)
      case(2)
         ! sphr = matmul(dtrafo, matmul(cart, transpose(dtrafo)))
         sphr(1, 1) = cart(3, 3) &
            & - 0.5_wp * (cart(3, 1) + cart(3, 2) + cart(1, 3) + cart(2, 3)) &
            & + 0.25_wp * (cart(1, 1) + cart(1, 2) + cart(2, 1) + cart(2, 2))
         sphr([2, 3, 5], 1) = s3 * cart([5, 6, 4], 3) &
            & - s3_4 * (cart([5, 6, 4], 1) + cart([5, 6, 4], 2))
         sphr(4, 1) = s3_4 * (cart(1, 3) - cart(2, 3)) &
            & - s3 * 0.25_wp * (cart(1, 1) - cart(2, 1) + cart(1, 2) - cart(2, 2))
         sphr(1, 2) = s3 * cart(3, 5) - s3_4 * (cart(1, 5) + cart(2, 5))
         sphr([2, 3, 5], 2) = 3 * cart([5, 6, 4], 5)
         sphr(4, 2) = 1.5_wp * (cart(1, 5) - cart(2, 5))
         sphr(1, 3) = s3 * cart(3, 6) - s3_4 * (cart(1, 6) + cart(2, 6))
         sphr([2, 3, 5], 3) = 3 * cart([5, 6, 4], 6)
         sphr(4, 3) = 1.5_wp * (cart(1, 6) - cart(2, 6))
         sphr(1, 4) = s3_4 * (cart(3, 1) - cart(3, 2)) &
            & - s3 * 0.25_wp * (cart(1, 1) - cart(1, 2) + cart(2, 1) - cart(2, 2))
         sphr([2, 3, 5], 4) = 1.5_wp * (cart([5, 6, 4], 1) - cart([5, 6, 4], 2))
         sphr(4, 4) = 0.75_wp * (cart(1, 1) - cart(2, 1) - cart(1, 2) + cart(2, 2))
         sphr(1, 5) = s3 * cart(3, 4) - s3_4 * (cart(1, 4) + cart(2, 4))
         sphr([2, 3, 5], 5) = 3 * cart([5, 6, 4], 4)
         sphr(4, 5) = 1.5_wp * (cart(1, 4) - cart(2, 4))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(dtrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(dtrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(3)
      select case(lj)
      case(0, 1)
         sphr = matmul(cart, transpose(ftrafo))
      case(2)
         sphr = matmul(dtrafo, matmul(cart, transpose(ftrafo)))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(ftrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(ftrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(4)
      select case(lj)
      case(0, 1)
         sphr = matmul(cart, transpose(gtrafo))
      case(2)
         sphr = matmul(dtrafo, matmul(cart, transpose(gtrafo)))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(gtrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(gtrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case default
      error stop "[Fatal] Moments higher than g are not supported"
   end select

end subroutine transform0

*/
/* TODO: priority low. Unforutnately this doesn't work */
/* Somehow, nvcc + template metaprogramming + a separate compilation file, gives an 'undefined reference to ' error. */
/* I will come back to this. For now the solution is to just copy the funcs to multipole.cu. */
// template <int A, int B>
// __device__
// void transform0(
//   const int lj, 
//   const int li, 
//   const float (&cart)[A][A], 
//         float (&sphr)[B][B]
// );

// template <int A, int B, int C>
// __device__
// void transform1(
//   const int lj, 
//   const int li, 
//   const float (&cart)[A][A][C], 
//         float (&sphr)[B][B][C]
// );
#endif