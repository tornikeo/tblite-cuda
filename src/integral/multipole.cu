#include "multipole.h"
#include "trafo.h"
#include <cassert>

template <int A, int B>
__device__
void transform0(
  const int lj, 
  const int li, 
  const float (&cart)[A][A], 
  float (&sphr)[B][B]
) {
  if (li == 0 || li == 1) {
    if (lj == 0 || lj == 1) {
      // Direct copy
      for (int i = 0; i < B; i++) {
        for (int j = 0; j < B; j++) {
          sphr[i][j] = cart[i][j];
        }
      }
    } else if (lj == 2) {
      // Transform for d-trafo
      for (int j = 0; j < B; j++) {
        sphr[0][j] = cart[2][j] - 0.5f * (cart[0][j] + cart[1][j]);
        sphr[1][j] = sqrtf(3.0f) * cart[4][j];
        sphr[2][j] = sqrtf(3.0f) * cart[5][j];
        sphr[3][j] = sqrtf(3.0f / 4.0f) * (cart[0][j] - cart[1][j]);
        sphr[4][j] = sqrtf(3.0f) * cart[3][j];
      }
    } else {
      // Unsupported moments higher than g
      assert(false && "[Fatal] Moments higher than g are not supported");
    }
  } else if (li == 2) {
    if (lj == 0 || lj == 1) {
      // Transform for d-trafo (transpose)
      for (int i = 0; i < B; i++) {
        sphr[i][0] = cart[i][2] - 0.5f * (cart[i][0] + cart[i][1]);
        sphr[i][1] = sqrtf(3.0f) * cart[i][4];
        sphr[i][2] = sqrtf(3.0f) * cart[i][5];
        sphr[i][3] = sqrtf(3.0f / 4.0f) * (cart[i][0] - cart[i][1]);
        sphr[i][4] = sqrtf(3.0f) * cart[i][3];
      }
    } else if (lj == 2) {
      // Transform for d-trafo (full)
      for (int i = 0; i < B; i++) {
        for (int j = 0; j < B; j++) {
          sphr[i][j] = cart[2][2] 
            - 0.5f * (cart[2][0] + cart[2][1] + cart[0][2] + cart[1][2]) 
            + 0.25f * (cart[0][0] + cart[0][1] + cart[1][0] + cart[1][1]);
        }
      }
    } else {
      // Unsupported moments higher than g
      assert(false && "[Fatal] Moments higher than g are not supported");
    }
  } else {
    // Unsupported moments higher than g
    assert(false && "[Fatal] Moments higher than g are not supported");
  }
}


template <int A, int B, int C>
__device__
void transform1(
  const int lj, 
  const int li, 
  const float (&cart)[A][A][C], /* TODO: priority medium. use better axes here, C should be first. */
        float (&sphr)[B][B][C]
) {
  for (int k = 0; k < C; ++k) {
    float tmpmat[B][B] = {0};

    for (int i = 0; i < B; ++i) {
      for (int j = 0; j < B; ++j) {
        tmpmat[i][j] = cart[k][i][j];
      }
    }

    transform0(lj, li, tmpmat, tmpmat);

    for (int i = 0; i < B; ++i) {
      for (int j = 0; j < B; ++j) {
        sphr[k][i][j] = tmpmat[i][j];
      }
    }
  }
}

template <int A, int B, int C>
__device__
void transform2(
  const int lj, 
  const int li, 
  const float (&cart)[A][A][C][C], /* TODO: priority medium. use better axes here, C should be first. */
        float (&sphr)[B][B][C][C]
) {
  for (int l = 0; l < A; ++l) {
    for (int k = 0; k < A; ++k) {
      float tmp_cart[C][C] = {0};
      float tmp_sphr[C][C] = {0};

      // Extract a 2D slice from the 4D cart array
      for (int i = 0; i < C; ++i) {
        for (int j = 0; j < C; ++j) {
          tmp_cart[i][j] = cart[k][l][i][j];
        }
      }

      // Perform the transform on the 2D slice
      transform0(lj, li, tmp_cart, tmp_sphr);

      // Store the transformed 2D slice back into the 4D sphr array
      for (int i = 0; i < C; ++i) {
        for (int j = 0; j < C; ++j) {
          sphr[k][l][i][j] = tmp_sphr[i][j];
        }
      }
    }
  }
}

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
__device__
inline float overlap_1d(const int moment, const int alpha)
{
  constexpr float dfactorial[] = {1.0f, 1.0f, 3.0f, 15.0f, 105.0f, 945.0f, 10395.0f, 135135.0f};
  if (moment % 2 == 0) {
    return powf(0.5f / alpha, moment / 2) * dfactorial[moment / 2];
  } else {
    return 0.0f;
  }
}

/*
pure subroutine horizontal_shift(ae, l, cfs)
   integer, intent(in) :: l
   real(wp), intent(in) :: ae
   real(wp), intent(inout) :: cfs(*)
   select case(l)
   case(0) ! s
      continue
   case(1) ! p
      cfs(1)=cfs(1)+ae*cfs(2)
   case(2) ! d
      cfs(1)=cfs(1)+ae*ae*cfs(3)
      cfs(2)=cfs(2)+ 2*ae*cfs(3)
   case(3) ! f
      cfs(1)=cfs(1)+ae*ae*ae*cfs(4)
      cfs(2)=cfs(2)+ 3*ae*ae*cfs(4)
      cfs(3)=cfs(3)+ 3*ae*cfs(4)
   case(4) ! g
      cfs(1)=cfs(1)+ae*ae*ae*ae*cfs(5)
      cfs(2)=cfs(2)+ 4*ae*ae*ae*cfs(5)
      cfs(3)=cfs(3)+ 6*ae*ae*cfs(5)
      cfs(4)=cfs(4)+ 4*ae*cfs(5)
   end select
end subroutine horizontal_shift
*/

__device__
void horizontal_shift(
  const float &ae, 
  const int &l, 
  float (&cfs)[])
{
  switch (l) {
    case 0: // s
      break;
    case 1: // p
      cfs[0] += ae * cfs[1];
      break;
    case 2: // d
      cfs[0] += ae * ae * cfs[2];
      cfs[1] += 2 * ae * cfs[2];
      break;
    case 3: // f
      cfs[0] += ae * ae * ae * cfs[3];
      cfs[1] += 3 * ae * ae * cfs[3];
      cfs[2] += 3 * ae * cfs[3];
      break;
    case 4: // g
      cfs[0] += ae * ae * ae * ae * cfs[4];
      cfs[1] += 4 * ae * ae * ae * cfs[4];
      cfs[2] += 6 * ae * ae * cfs[4];
      cfs[3] += 4 * ae * cfs[4];
      break;
  }
}

/*

pure subroutine form_product(a, b, la, lb, d)
   integer, intent(in) :: la, lb
   real(wp), intent(in) :: a(*), b(*)
   real(wp), intent(inout) :: d(*)
   if(la.ge.4.or.lb.ge.4) goto 40
   if(la.ge.3.or.lb.ge.3) goto 30
   if(la.ge.2.or.lb.ge.2) goto 20
   ! <s|s> = <s>
   d(1)=a(1)*b(1)
   if(la.eq.0.and.lb.eq.0) return
   ! <s|p> = <s|*(|s>+|p>)
   !       = <s> + <p>
   d(2)=a(1)*b(2)+a(2)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|p> = (<s|+<p|)*(|s>+|p>)
   !       = <s> + <p> + <d>
   d(3)=a(2)*b(2)
   return
20 continue
   ! <s|d> = <s|*(|s>+|p>+|d>)
   !       = <s> + <p> + <d>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f>
   d(3)=d(3)+a(2)*b(2)
   d(4)=a(2)*b(3)+a(3)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(5)=a(3)*b(3)
   return
30 continue
   ! <s|f> = <s|*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=a(2)*b(4)+a(4)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(5)=d(5)+a(3)*b(3)
   d(6)=a(3)*b(4)+a(4)*b(3)
   if(la.le.2.or.lb.le.2) return
   ! <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(7)=a(4)*b(4)
   return
40 continue
   ! <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   d(5)=a(1)*b(5)+a(5)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=d(5)+a(2)*b(4)+a(4)*b(2)
   d(6)=a(2)*b(5)+a(5)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(5)=d(5)+a(3)*b(3)
   d(6)=d(5)+a(3)*b(4)+a(4)*b(3)
   d(7)=a(3)*b(5)+a(5)*b(3)
   if(la.le.2.or.lb.le.2) return
   ! <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
   d(7)=d(7)+a(4)*b(4)
   d(8)=a(4)*b(5)+a(5)*b(4)
   if(la.le.3.or.lb.le.3) return
   ! <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
   d(9)=a(5)*b(5)

end subroutine form_product

*/




__device__
/* This function is repulsive (to me) */
inline void form_product(
  const float (&a)[],
  const float (&b)[],
  const int la, 
  const int lb,
  float (&d)[]
) {
  if (la >= 4 || lb >= 4) {
    // <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
    //       = <s> + <p> + <d> + <f> + <g>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    d[3] = a[0] * b[3] + a[3] * b[0];
    d[4] = a[0] * b[4] + a[4] * b[0];
    if (la == 0 || lb == 0) return;

    // <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
    //       = <s> + <p> + <d> + <f> + <g> + <h>
    d[2] += a[1] * b[1];
    d[3] += a[1] * b[2] + a[2] * b[1];
    d[4] += a[1] * b[3] + a[3] * b[1];
    d[5] = a[1] * b[4] + a[4] * b[1];
    if (la <= 1 || lb <= 1) return;

    // <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
    //       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
    d[4] += a[2] * b[2];
    d[5] += a[2] * b[3] + a[3] * b[2];
    d[6] = a[2] * b[4] + a[4] * b[2];
    if (la <= 2 || lb <= 2) return;

    // <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
    //       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
    d[6] += a[3] * b[3];
    d[7] = a[3] * b[4] + a[4] * b[3];
    if (la <= 3 || lb <= 3) return;

    // <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
    //       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
    d[8] = a[4] * b[4];
    return;
  }

  if (la >= 3 || lb >= 3) {
    // <s|f> = <s|*(|s>+|p>+|d>+|f>)
    //       = <s> + <p> + <d> + <f>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    d[3] = a[0] * b[3] + a[3] * b[0];
    if (la == 0 || lb == 0) return;

    // <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
    //       = <s> + <p> + <d> + <f> + <g>
    d[2] += a[1] * b[1];
    d[3] += a[1] * b[2] + a[2] * b[1];
    d[4] = a[1] * b[3] + a[3] * b[1];
    if (la <= 1 || lb <= 1) return;

    // <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
    //       = <s> + <p> + <d> + <f> + <g> + <h>
    d[4] += a[2] * b[2];
    d[5] = a[2] * b[3] + a[3] * b[2];
    if (la <= 2 || lb <= 2) return;

    // <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
    //       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
    d[6] = a[3] * b[3];
    return;
  }

  if (la >= 2 || lb >= 2) {
    // <s|d> = <s|*(|s>+|p>+|d>)
    //       = <s> + <p> + <d>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    if (la == 0 || lb == 0) return;

    // <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
    //       = <s> + <p> + <d> + <f>
    d[2] += a[1] * b[1];
    d[3] = a[1] * b[2] + a[2] * b[1];
    if (la <= 1 || lb <= 1) return;

    // <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
    //       = <s> + <p> + <d> + <f> + <g>
    d[4] = a[2] * b[2];
    return;
  }

  // <s|s> = <s>
  d[0] = a[0] * b[0];
  if (la == 0 && lb == 0) return;

  // <s|p> = <s|*(|s>+|p>)
  //       = <s> + <p>
  d[1] = a[0] * b[1] + a[1] * b[0];
  if (la == 0 || lb == 0) return;

  // <p|p> = (<s|+<p|)*(|s>+|p>)
  //       = <s> + <p> + <d>
  d[2] = a[1] * b[1];
}

/*
pure subroutine multipole_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, d3d, q3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: d3d(3)
   real(wp), intent(out) :: q3d(6)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3, 3)

   v1d(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         v1d(k, 1) = v1d(k, 1) + s1d(l) * vv(l)
         v1d(k, 2) = v1d(k, 2) + (s1d(l+1) + rpi(k)*s1d(l)) * vv(l)
         v1d(k, 3) = v1d(k, 3) + (s1d(l+2) + 2*rpi(k)*s1d(l+1) + rpi(k)*rpi(k)*s1d(l)) * vv(l)
      end do
   end do

   s3d = v1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   d3d(1) = v1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   d3d(2) = v1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   d3d(3) = v1d(1, 1) * v1d(2, 1) * v1d(3, 2)
   q3d(1) = v1d(1, 3) * v1d(2, 1) * v1d(3, 1)
   q3d(2) = v1d(1, 2) * v1d(2, 2) * v1d(3, 1)
   q3d(3) = v1d(1, 1) * v1d(2, 3) * v1d(3, 1)
   q3d(4) = v1d(1, 2) * v1d(2, 1) * v1d(3, 2)
   q3d(5) = v1d(1, 1) * v1d(2, 2) * v1d(3, 2)
   q3d(6) = v1d(1, 1) * v1d(2, 1) * v1d(3, 3)

end subroutine multipole_3d
*/

__device__ 
void multipole_3d(
  const float (&rpi)[3],
  const float (&rpj)[3],
  const float ai,
  const float aj,
  const int (&li)[3],
  const int (&lj)[3],
  const float (&s1d)[MAXL2],
  float &s3d,
  float (&d3d)[3],
  float (&q3d)[6]
) {
  float vi[MAXL] = {0}, vj[MAXL] = {0}, vv[MAXL2] = {0}, v1d[3][3] = {0};

  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < MAXL; ++i) {
      vi[i] = 0.0f;
      vj[i] = 0.0f;
    }
    
    vi[li[k]] = 1.0f;
    vj[lj[k]] = 1.0f;
    
    horizontal_shift(rpi[k], li[k], vi);
    horizontal_shift(rpj[k], lj[k], vj);
    form_product(vi, vj, li[k], lj[k], vv);

    for (int l = 0; l <= li[k] + lj[k]; ++l) {
      v1d[k][0] += s1d[l] * vv[l];
      v1d[k][1] += (s1d[l + 1] + rpi[k] * s1d[l]) * vv[l];
      v1d[k][2] += (s1d[l + 2] + 2 * rpi[k] * s1d[l + 1] + rpi[k] * rpi[k] * s1d[l]) * vv[l];
    }
  }

  s3d = v1d[0][0] * v1d[1][0] * v1d[2][0];
  d3d[0] = v1d[0][1] * v1d[1][0] * v1d[2][0];
  d3d[1] = v1d[0][0] * v1d[1][1] * v1d[2][0];
  d3d[2] = v1d[0][0] * v1d[1][0] * v1d[2][1];
  q3d[0] = v1d[0][2] * v1d[1][0] * v1d[2][0];
  q3d[1] = v1d[0][1] * v1d[1][1] * v1d[2][0];
  q3d[2] = v1d[0][0] * v1d[1][2] * v1d[2][0];
  q3d[3] = v1d[0][1] * v1d[1][0] * v1d[2][1];
  q3d[4] = v1d[0][0] * v1d[1][1] * v1d[2][1];
  q3d[5] = v1d[0][0] * v1d[1][0] * v1d[2][2];
}




__device__
void multipole_cgto(
  const cgto_type &cgtoi,
  const cgto_type &cgtoj,
  const float r2,
  const float (&vec)[3], /* TODO: Priority low, this could be a builtin float3 */
  const float intcut,
  float (&overlap)[msao[MAXL]][msao[MAXL]],
  float (&dpint  )[msao[MAXL]][msao[MAXL]][3],
  float (&qpint  )[msao[MAXL]][msao[MAXL]][6]
)
{
  /* TODO: This feels wrong, and is it?
     NOTE: I tried moving this outside as an operator[] arr, see msao
     doesn't work nicely, because 2D arr indexing needs more code. This 
     is preferred for now.
  */
  constexpr int lx[mlao(MAXL) + lmap(MAXL)][3] = {
    {0,0,0},
    {0,1,0},
    {0,0,1},
    {1,0,0},
    {2,0,0},
    {0,2,0},
    {0,0,2},
    {1,1,0},
    {1,0,1},
    {0,1,1}
    // {3,0,0},
    // {0,3,0},
    // {0,0,3},
    // {2,1,0},
    // {2,0,1},
    // {1,2,0},
    // {0,2,1},
    // {1,0,2},
    // {0,1,2},
    // {1,1,1},
    // {4,0,0},
    // {0,4,0},
    // {0,0,4},
    // {3,1,0}
  };
  /*
  integer :: ip, jp, mli, mlj, l
  real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6), pre, tr
  real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
  real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
  real(wp) :: q3d(6, mlao(cgtoj%ang), mlao(cgtoi%ang))

  s3d(:, :) = 0.0_wp
  d3d(:, :, :) = 0.0_wp
  q3d(:, :, :) = 0.0_wp
  */
  float ip, jp, mli, mlj, l;
  float eab, oab, est, s1d[MAXL2], rpi[3], rpj[3], cc, val, dip[3], quad[6], pre, tr;
  float s3d[mlao(MAXL)][mlao(MAXL)] = {0};
  float d3d[mlao(MAXL)][mlao(MAXL)][3] = {0};
  float q3d[mlao(MAXL)][mlao(MAXL)][6] = {0};
  /*
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
  */

  for (int ip = 0; ip < cgtoi.nprim; ++ip) {
    for (int jp = 0; jp < cgtoj.nprim; ++jp) {
      eab = cgtoi.alpha[ip] + cgtoj.alpha[jp];
      oab = 1.0f / eab;
      est = cgtoi.alpha[ip] * cgtoj.alpha[jp] * r2 * oab;
      if (est > intcut) continue;
      pre = expf(-est) * sqrtf(M_PI * M_PI * M_PI) * powf(oab, 1.5f);
      for (int i = 0; i < 3; ++i) {
        rpi[i] = -vec[i] * cgtoj.alpha[jp] * oab;
        rpj[i] = vec[i] * cgtoi.alpha[ip] * oab;
      }
      for (int l = 0; l <= cgtoi.ang + cgtoj.ang + 2; ++l) {
        s1d[l] = overlap_1d(l, eab);
      }
      cc = cgtoi.coeff[ip] * cgtoj.coeff[jp] * pre;
      for (int mli = 0; mli < mlao(cgtoi.ang); ++mli) {
        for (int mlj = 0; mlj < mlao(cgtoj.ang); ++mlj) {
          float val, dip[3], quad[6];
          multipole_3d(
            rpj, rpi, cgtoj.alpha[jp], cgtoi.alpha[ip],
            lx[mlj + lmap(cgtoj.ang)], lx[mli + lmap(cgtoi.ang)],
            s1d, val, dip, quad
          );
          s3d[mlj][mli] += cc * val;
          for (int i = 0; i < 3; ++i) {
            d3d[i][mlj][mli] += cc * dip[i];
          }
          for (int i = 0; i < 6; ++i) {
            q3d[i][mlj][mli] += cc * quad[i];
          }
        }
      }
    }
  }

  transform0(cgtoj.ang, cgtoi.ang, s3d, overlap);
  transform1(cgtoj.ang, cgtoi.ang, d3d, dpint);
  transform1(cgtoj.ang, cgtoi.ang, q3d, qpint);

  /*   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do*/
  for (int mli = 0; mli < msao[cgtoi.ang]; mli++)
  {
    for (int mlj = 0; mlj < msao[cgtoj.ang]; mlj++)
    {
      float tr = 0.5f * (qpint[mlj][mli][0] + qpint[mlj][mli][2] + qpint[mlj][mli][5]);
      qpint[mlj][mli][0] = 1.5f * qpint[mlj][mli][0] - tr;
      qpint[mlj][mli][1] = 1.5f * qpint[mlj][mli][1];
      qpint[mlj][mli][2] = 1.5f * qpint[mlj][mli][2] - tr;
      qpint[mlj][mli][3] = 1.5f * qpint[mlj][mli][3];
      qpint[mlj][mli][4] = 1.5f * qpint[mlj][mli][4];
      qpint[mlj][mli][5] = 1.5f * qpint[mlj][mli][5] - tr;
    }
  }
}