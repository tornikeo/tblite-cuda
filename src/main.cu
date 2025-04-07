#include "main.h"
#include <iostream>
#include <assert.h>
// pure subroutine form_product(a, b, la, lb, d)
//    integer, intent(in) :: la, lb
//    real(wp), intent(in) :: a(*), b(*)
//    real(wp), intent(inout) :: d(*)
//    if(la.ge.4.or.lb.ge.4) goto 40
//    if(la.ge.3.or.lb.ge.3) goto 30
//    if(la.ge.2.or.lb.ge.2) goto 20
//    ! <s|s> = <s>
//    d(1)=a(1)*b(1)
//    if(la.eq.0.and.lb.eq.0) return
//    ! <s|p> = <s|*(|s>+|p>)
//    !       = <s> + <p>
//    d(2)=a(1)*b(2)+a(2)*b(1)
//    if(la.eq.0.or.lb.eq.0) return
//    ! <p|p> = (<s|+<p|)*(|s>+|p>)
//    !       = <s> + <p> + <d>
//    d(3)=a(2)*b(2)
//    return
// 20 continue
//    ! <s|d> = <s|*(|s>+|p>+|d>)
//    !       = <s> + <p> + <d>
//    d(1)=a(1)*b(1)
//    d(2)=a(1)*b(2)+a(2)*b(1)
//    d(3)=a(1)*b(3)+a(3)*b(1)
//    if(la.eq.0.or.lb.eq.0) return
//    ! <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
//    !       = <s> + <p> + <d> + <f>
//    d(3)=d(3)+a(2)*b(2)
//    d(4)=a(2)*b(3)+a(3)*b(2)
//    if(la.le.1.or.lb.le.1) return
//    ! <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
//    !       = <s> + <p> + <d> + <f> + <g>
//    d(5)=a(3)*b(3)
//    return
// 30 continue
//    ! <s|f> = <s|*(|s>+|p>+|d>+|f>)
//    !       = <s> + <p> + <d> + <f>
//    d(1)=a(1)*b(1)
//    d(2)=a(1)*b(2)+a(2)*b(1)
//    d(3)=a(1)*b(3)+a(3)*b(1)
//    d(4)=a(1)*b(4)+a(4)*b(1)
//    if(la.eq.0.or.lb.eq.0) return
//    ! <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
//    !       = <s> + <p> + <d> + <f> + <g>
//    d(3)=d(3)+a(2)*b(2)
//    d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
//    d(5)=a(2)*b(4)+a(4)*b(2)
//    if(la.le.1.or.lb.le.1) return
//    ! <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
//    !       = <s> + <p> + <d> + <f> + <g> + <h>
//    d(5)=d(5)+a(3)*b(3)
//    d(6)=a(3)*b(4)+a(4)*b(3)
//    if(la.le.2.or.lb.le.2) return
//    ! <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
//    !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
//    d(7)=a(4)*b(4)
//    return
// 40 continue
//    ! <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
//    !       = <s> + <p> + <d> + <f> + <g>
//    d(1)=a(1)*b(1)
//    d(2)=a(1)*b(2)+a(2)*b(1)
//    d(3)=a(1)*b(3)+a(3)*b(1)
//    d(4)=a(1)*b(4)+a(4)*b(1)
//    d(5)=a(1)*b(5)+a(5)*b(1)
//    if(la.eq.0.or.lb.eq.0) return
//    ! <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
//    !       = <s> + <p> + <d> + <f> + <g> + <h>
//    d(3)=d(3)+a(2)*b(2)
//    d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
//    d(5)=d(5)+a(2)*b(4)+a(4)*b(2)
//    d(6)=a(2)*b(5)+a(5)*b(2)
//    if(la.le.1.or.lb.le.1) return
//    ! <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
//    !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
//    d(5)=d(5)+a(3)*b(3)
//    d(6)=d(5)+a(3)*b(4)+a(4)*b(3)
//    d(7)=a(3)*b(5)+a(5)*b(3)
//    if(la.le.2.or.lb.le.2) return
//    ! <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
//    !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
//    d(7)=d(7)+a(4)*b(4)
//    d(8)=a(4)*b(5)+a(5)*b(4)
//    if(la.le.3.or.lb.le.3) return
//    ! <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
//    !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
//    d(9)=a(5)*b(5)

// end subroutine form_product

__device__ void form_product(const double *a, const double *b, int la, int lb, double *d)
{
    if (la >= 4 || lb >= 4)
        goto level_40;
    if (la >= 3 || lb >= 3)
        goto level_30;
    if (la >= 2 || lb >= 2)
        goto level_20;

    // <s|s> = <s>
    d[0] = a[0] * b[0];
    if (la == 0 && lb == 0)
        return;

    // <s|p> = <s|*(|s>+|p>) = <s> + <p>
    d[1] = a[0] * b[1] + a[1] * b[0];
    if (la == 0 || lb == 0)
        return;

    // <p|p> = (<s|+<p|)*(|s>+|p>) = <s> + <p> + <d>
    d[2] = a[1] * b[1];
    return;

level_20:
    // <s|d> = <s|*(|s>+|p>+|d>) = <s> + <p> + <d>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    if (la == 0 || lb == 0)
        return;

    // <p|d> = (<s|+<p|)*(|s>+|p>+|d>) = <s> + <p> + <d> + <f>
    d[2] += a[1] * b[1];
    d[3] = a[1] * b[2] + a[2] * b[1];
    if (la <= 1 || lb <= 1)
        return;

    // <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>) = <s> + <p> + <d> + <f> + <g>
    d[4] = a[2] * b[2];
    return;

level_30:
    // <s|f> = <s|*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    d[3] = a[0] * b[3] + a[3] * b[0];
    if (la == 0 || lb == 0)
        return;

    // <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g>
    d[2] += a[1] * b[1];
    d[3] += a[1] * b[2] + a[2] * b[1];
    d[4] = a[1] * b[3] + a[3] * b[1];
    if (la <= 1 || lb <= 1)
        return;

    // <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g> + <h>
    d[4] += a[2] * b[2];
    d[5] = a[2] * b[3] + a[3] * b[2];
    if (la <= 2 || lb <= 2)
        return;

    // <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
    d[6] = a[3] * b[3];
    return;

level_40:
    // <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    d[3] = a[0] * b[3] + a[3] * b[0];
    d[4] = a[0] * b[4] + a[4] * b[0];
    if (la == 0 || lb == 0)
        return;

    // <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h>
    d[2] += a[1] * b[1];
    d[3] += a[1] * b[2] + a[2] * b[1];
    d[4] += a[1] * b[3] + a[3] * b[1];
    d[5] = a[1] * b[4] + a[4] * b[1];
    if (la <= 1 || lb <= 1)
        return;

    // <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
    d[4] += a[2] * b[2];
    d[5] += a[2] * b[3] + a[3] * b[2];
    d[6] = a[2] * b[4] + a[4] * b[2];
    if (la <= 2 || lb <= 2)
        return;

    // <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
    d[6] += a[3] * b[3];
    d[7] = a[3] * b[4] + a[4] * b[3];
    if (la <= 3 || lb <= 3)
        return;

    // <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
    d[8] = a[4] * b[4];
}

// pure subroutine horizontal_shift(ae, l, cfs)
//    integer, intent(in) :: l
//    real(wp), intent(in) :: ae
//    real(wp), intent(inout) :: cfs(*)
//    select case(l)
//    case(0) ! s
//       continue
//    case(1) ! p
//       cfs(1)=cfs(1)+ae*cfs(2)
//    case(2) ! d
//       cfs(1)=cfs(1)+ae*ae*cfs(3)
//       cfs(2)=cfs(2)+ 2*ae*cfs(3)
//    case(3) ! f
//       cfs(1)=cfs(1)+ae*ae*ae*cfs(4)
//       cfs(2)=cfs(2)+ 3*ae*ae*cfs(4)
//       cfs(3)=cfs(3)+ 3*ae*cfs(4)
//    case(4) ! g
//       cfs(1)=cfs(1)+ae*ae*ae*ae*cfs(5)
//       cfs(2)=cfs(2)+ 4*ae*ae*ae*cfs(5)
//       cfs(3)=cfs(3)+ 6*ae*ae*cfs(5)
//       cfs(4)=cfs(4)+ 4*ae*cfs(5)
//    end select
// end subroutine horizontal_shift

__device__ void horizontal_shift(double ae, int l, double *cfs)
{
    switch (l)
    {
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
    default:
        break;
    }
}

// ddi is [3, 3]
// dqi is [3, 6]
// ddj is [3, 3]
// dqj is [3, 6]
// pure subroutine shift_operator(vec, s, di, qi, ds, ddi, dqi, ddj, dqj)
//    real(wp),intent(in) :: vec(:)
//    real(wp),intent(in) :: s
//    real(wp),intent(in) :: di(:)
//    real(wp),intent(in) :: qi(:)
//    real(wp),intent(in) :: ds(:)
//    real(wp),intent(in) :: ddi(:, :)
//    real(wp),intent(in) :: dqi(:, :)
//    real(wp),intent(out) :: ddj(:, :)
//    real(wp),intent(out) :: dqj(:, :)

//    ddj(:, 1) = ddi(:, 1) - vec(1)*ds
//    ddj(:, 2) = ddi(:, 2) - vec(2)*ds
//    ddj(:, 3) = ddi(:, 3) - vec(3)*ds
//    ddj(1, 1) = ddj(1, 1) - s
//    ddj(2, 2) = ddj(2, 2) - s
//    ddj(3, 3) = ddj(3, 3) - s

//    dqj(:, 1) = dqi(:, 1) - 2*vec(1)*ddi(:, 1) + vec(1)**2*ds
//    dqj(:, 3) = dqi(:, 3) - 2*vec(2)*ddi(:, 2) + vec(2)**2*ds
//    dqj(:, 6) = dqi(:, 6) - 2*vec(3)*ddi(:, 3) + vec(3)**2*ds
//    dqj(:, 2) = dqi(:, 2) - vec(1)*ddi(:, 2) - vec(2)*ddi(:, 1) + vec(1)*vec(2)*ds
//    dqj(:, 4) = dqi(:, 4) - vec(1)*ddi(:, 3) - vec(3)*ddi(:, 1) + vec(1)*vec(3)*ds
//    dqj(:, 5) = dqi(:, 5) - vec(2)*ddi(:, 3) - vec(3)*ddi(:, 2) + vec(2)*vec(3)*ds
//    dqj(1, 1) = dqj(1, 1) - 2*di(1) + 2*vec(1)*s
//    dqj(2, 3) = dqj(2, 3) - 2*di(2) + 2*vec(2)*s
//    dqj(3, 6) = dqj(3, 6) - 2*di(3) + 2*vec(3)*s
//    dqj(1, 2) = dqj(1, 2) - di(2) + vec(2)*s
//    dqj(2, 2) = dqj(2, 2) - di(1) + vec(1)*s
//    dqj(1, 4) = dqj(1, 4) - di(3) + vec(3)*s
//    dqj(3, 4) = dqj(3, 4) - di(1) + vec(1)*s
//    dqj(2, 5) = dqj(2, 5) - di(3) + vec(3)*s
//    dqj(3, 5) = dqj(3, 5) - di(2) + vec(2)*s

// end subroutine shift_operator

__device__ void shift_operator(const double *vec, double s, const double *di, const double *qi, const double *ds,
                               const double ddi[3][3], const double dqi[3][6], double ddj[3][3], double dqj[3][6])
{
    // Update ddj
    for (int i = 0; i < 3; i++)
    {
        ddj[i][0] = ddi[i][0] - vec[0] * ds[i];
        ddj[i][1] = ddi[i][1] - vec[1] * ds[i];
        ddj[i][2] = ddi[i][2] - vec[2] * ds[i];
    }
    ddj[0][0] -= s;
    ddj[1][1] -= s;
    ddj[2][2] -= s;

    // Update dqj
    for (int i = 0; i < 3; i++)
    {
        dqj[i][0] = dqi[i][0] - 2 * vec[0] * ddi[i][0] + vec[0] * vec[0] * ds[i];
        dqj[i][2] = dqi[i][2] - 2 * vec[1] * ddi[i][1] + vec[1] * vec[1] * ds[i];
        dqj[i][5] = dqi[i][5] - 2 * vec[2] * ddi[i][2] + vec[2] * vec[2] * ds[i];
        dqj[i][1] = dqi[i][1] - vec[0] * ddi[i][1] - vec[1] * ddi[i][0] + vec[0] * vec[1] * ds[i];
        dqj[i][3] = dqi[i][3] - vec[0] * ddi[i][2] - vec[2] * ddi[i][0] + vec[0] * vec[2] * ds[i];
        dqj[i][4] = dqi[i][4] - vec[1] * ddi[i][2] - vec[2] * ddi[i][1] + vec[1] * vec[2] * ds[i];
    }
    dqj[0][0] -= 2 * di[0] - 2 * vec[0] * s;
    dqj[1][2] -= 2 * di[1] - 2 * vec[1] * s;
    dqj[2][5] -= 2 * di[2] - 2 * vec[2] * s;
    dqj[0][1] -= di[1] - vec[1] * s;
    dqj[1][1] -= di[0] - vec[0] * s;
    dqj[0][3] -= di[2] - vec[2] * s;
    dqj[2][3] -= di[0] - vec[0] * s;
    dqj[1][4] -= di[2] - vec[2] * s;
    dqj[2][4] -= di[1] - vec[1] * s;
}


// pure subroutine multipole_grad_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, d3d, q3d, &
//     & ds3d, dd3d, dq3d)
//  real(wp), intent(in) :: rpi(3)
//  real(wp), intent(in) :: rpj(3)
//  real(wp), intent(in) :: ai
//  real(wp), intent(in) :: aj
//  integer, intent(in) :: li(3)
//  integer, intent(in) :: lj(3)
//  real(wp), intent(in) :: s1d(0:)
//  real(wp), intent(out) :: s3d
//  real(wp), intent(out) :: d3d(3)
//  real(wp), intent(out) :: q3d(6)
//  real(wp), intent(out) :: ds3d(3)
//  real(wp), intent(out) :: dd3d(3, 3)
//  real(wp), intent(out) :: dq3d(3, 6)

//  integer :: k, l
//  real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3, 3)
//  real(wp) :: gi(0:maxl), gg(0:maxl2), g1d(3, 3), rpc

//  v1d(:, :) = 0.0_wp
//  g1d(:, :) = 0.0_wp

//  do k = 1, 3
//     vv(:) = 0.0_wp
//     gg(:) = 0.0_wp
//     vi(:) = 0.0_wp
//     vj(:) = 0.0_wp
//     gi(:) = 0.0_wp
//     rpc = rpj(k)

//     vi(li(k)) = 1.0_wp
//     vj(lj(k)) = 1.0_wp
//     gi(li(k)+1) = 2*ai
//     if (li(k) > 0) gi(li(k)-1) = -li(k)

//     call horizontal_shift(rpi(k), li(k)-1, gi)
//     call horizontal_shift(rpi(k), li(k)+1, gi)
//     call horizontal_shift(rpi(k), li(k), vi)
//     call horizontal_shift(rpj(k), lj(k), vj)
//     call form_product(vi, vj, li(k), lj(k), vv)
//     call form_product(gi, vj, li(k)+1, lj(k), gg)
//     do l = 0, li(k) + lj(k) + 1
//        v1d(k, 1) = v1d(k, 1) + s1d(l) * vv(l)
//        v1d(k, 2) = v1d(k, 2) + (s1d(l+1) + rpc*s1d(l)) * vv(l)
//        v1d(k, 3) = v1d(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * vv(l)
//        g1d(k, 1) = g1d(k, 1) + s1d(l) * gg(l)
//        g1d(k, 2) = g1d(k, 2) + (s1d(l+1) + rpc*s1d(l)) * gg(l)
//        g1d(k, 3) = g1d(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * gg(l)
//     end do
//  end do

//  s3d = v1d(1, 1) * v1d(2, 1) * v1d(3, 1)
//  d3d(1) = v1d(1, 2) * v1d(2, 1) * v1d(3, 1)
//  d3d(2) = v1d(1, 1) * v1d(2, 2) * v1d(3, 1)
//  d3d(3) = v1d(1, 1) * v1d(2, 1) * v1d(3, 2)
//  q3d(1) = v1d(1, 3) * v1d(2, 1) * v1d(3, 1)
//  q3d(2) = v1d(1, 2) * v1d(2, 2) * v1d(3, 1)
//  q3d(3) = v1d(1, 1) * v1d(2, 3) * v1d(3, 1)
//  q3d(4) = v1d(1, 2) * v1d(2, 1) * v1d(3, 2)
//  q3d(5) = v1d(1, 1) * v1d(2, 2) * v1d(3, 2)
//  q3d(6) = v1d(1, 1) * v1d(2, 1) * v1d(3, 3)

//  ds3d(1) = g1d(1, 1) * v1d(2, 1) * v1d(3, 1)
//  ds3d(2) = v1d(1, 1) * g1d(2, 1) * v1d(3, 1)
//  ds3d(3) = v1d(1, 1) * v1d(2, 1) * g1d(3, 1)
//  dd3d(1, 1) = g1d(1, 2) * v1d(2, 1) * v1d(3, 1)
//  dd3d(2, 1) = v1d(1, 2) * g1d(2, 1) * v1d(3, 1)
//  dd3d(3, 1) = v1d(1, 2) * v1d(2, 1) * g1d(3, 1)
//  dd3d(1, 2) = g1d(1, 1) * v1d(2, 2) * v1d(3, 1)
//  dd3d(2, 2) = v1d(1, 1) * g1d(2, 2) * v1d(3, 1)
//  dd3d(3, 2) = v1d(1, 1) * v1d(2, 2) * g1d(3, 1)
//  dd3d(1, 3) = g1d(1, 1) * v1d(2, 1) * v1d(3, 2)
//  dd3d(2, 3) = v1d(1, 1) * g1d(2, 1) * v1d(3, 2)
//  dd3d(3, 3) = v1d(1, 1) * v1d(2, 1) * g1d(3, 2)
//  dq3d(1, 1) = g1d(1, 3) * v1d(2, 1) * v1d(3, 1)
//  dq3d(2, 1) = v1d(1, 3) * g1d(2, 1) * v1d(3, 1)
//  dq3d(3, 1) = v1d(1, 3) * v1d(2, 1) * g1d(3, 1)
//  dq3d(1, 2) = g1d(1, 2) * v1d(2, 2) * v1d(3, 1)
//  dq3d(2, 2) = v1d(1, 2) * g1d(2, 2) * v1d(3, 1)
//  dq3d(3, 2) = v1d(1, 2) * v1d(2, 2) * g1d(3, 1)
//  dq3d(1, 3) = g1d(1, 1) * v1d(2, 3) * v1d(3, 1)
//  dq3d(2, 3) = v1d(1, 1) * g1d(2, 3) * v1d(3, 1)
//  dq3d(3, 3) = v1d(1, 1) * v1d(2, 3) * g1d(3, 1)
//  dq3d(1, 4) = g1d(1, 2) * v1d(2, 1) * v1d(3, 2)
//  dq3d(2, 4) = v1d(1, 2) * g1d(2, 1) * v1d(3, 2)
//  dq3d(3, 4) = v1d(1, 2) * v1d(2, 1) * g1d(3, 2)
//  dq3d(1, 5) = g1d(1, 1) * v1d(2, 2) * v1d(3, 2)
//  dq3d(2, 5) = v1d(1, 1) * g1d(2, 2) * v1d(3, 2)
//  dq3d(3, 5) = v1d(1, 1) * v1d(2, 2) * g1d(3, 2)
//  dq3d(1, 6) = g1d(1, 1) * v1d(2, 1) * v1d(3, 3)
//  dq3d(2, 6) = v1d(1, 1) * g1d(2, 1) * v1d(3, 3)
//  dq3d(3, 6) = v1d(1, 1) * v1d(2, 1) * g1d(3, 3)

// end subroutine multipole_grad_3d

__device__ 
void multipole_grad_3d(
    const double rpi[3], const double rpj[3],
    const double ai, const double aj, 
    const int li[3], const int lj[3], const double s1d[MAXL2],
    
    double &s3d, double d3d[3], double q3d[3], double ds3d[3], 
    double dd3d[3][3], double dq3d[3][6])
{
    double v1d[3][3] = {0.0};
    double g1d[3][3] = {0.0};
    double vi[MAXL + 1] = {0.0};
    double vj[MAXL + 1] = {0.0};
    double vv[MAXL2 + 1] = {0.0};
    double gi[MAXL + 1] = {0.0};
    double gg[MAXL2 + 1] = {0.0};

    for (int k = 0; k < 3; k++)
    {
        double rpc = rpj[k];
        for (int i = 0; i <= MAXL; i++)
        {
            vi[i] = 0.0;
            vj[i] = 0.0;
            gi[i] = 0.0;
        }
        for (int i = 0; i <= MAXL2; i++)
        {
            vv[i] = 0.0;
            gg[i] = 0.0;
        }

        vi[li[k]] = 1.0;
        vj[lj[k]] = 1.0;
        gi[li[k] + 1] = 2 * ai;
        if (li[k] > 0)
            gi[li[k] - 1] = -li[k];

        horizontal_shift(rpi[k], li[k] - 1, gi);
        horizontal_shift(rpi[k], li[k] + 1, gi);
        horizontal_shift(rpi[k], li[k], vi);
        horizontal_shift(rpj[k], lj[k], vj);
        form_product(vi, vj, li[k], lj[k], vv);
        form_product(gi, vj, li[k] + 1, lj[k], gg);

        for (int l = 0; l <= li[k] + lj[k] + 1; l++)
        {
            v1d[k][0] += s1d[l] * vv[l];
            v1d[k][1] += (s1d[l + 1] + rpc * s1d[l]) * vv[l];
            v1d[k][2] += (s1d[l + 2] + 2 * rpc * s1d[l + 1] + rpc * rpc * s1d[l]) * vv[l];
            g1d[k][0] += s1d[l] * gg[l];
            g1d[k][1] += (s1d[l + 1] + rpc * s1d[l]) * gg[l];
            g1d[k][2] += (s1d[l + 2] + 2 * rpc * s1d[l + 1] + rpc * rpc * s1d[l]) * gg[l];
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

    ds3d[0] = g1d[0][0] * v1d[1][0] * v1d[2][0];
    ds3d[1] = v1d[0][0] * g1d[1][0] * v1d[2][0];
    ds3d[2] = v1d[0][0] * v1d[1][0] * g1d[2][0];
    dd3d[0][0] = g1d[0][1] * v1d[1][0] * v1d[2][0];
    dd3d[1][0] = v1d[0][1] * g1d[1][0] * v1d[2][0];
    dd3d[2][0] = v1d[0][1] * v1d[1][0] * g1d[2][0];
    dd3d[0][1] = g1d[0][0] * v1d[1][1] * v1d[2][0];
    dd3d[1][1] = v1d[0][0] * g1d[1][1] * v1d[2][0];
    dd3d[2][1] = v1d[0][0] * v1d[1][1] * g1d[2][0];
    dd3d[0][2] = g1d[0][0] * v1d[1][0] * v1d[2][1];
    dd3d[1][2] = v1d[0][0] * g1d[1][0] * v1d[2][1];
    dd3d[2][2] = v1d[0][0] * v1d[1][0] * g1d[2][1];
    dq3d[0][0] = g1d[0][2] * v1d[1][0] * v1d[2][0];
    dq3d[1][0] = v1d[0][2] * g1d[1][0] * v1d[2][0];
    dq3d[2][0] = v1d[0][2] * v1d[1][0] * g1d[2][0];
    dq3d[0][1] = g1d[0][1] * v1d[1][1] * v1d[2][0];
    dq3d[1][1] = v1d[0][1] * g1d[1][1] * v1d[2][0];
    dq3d[2][1] = v1d[0][1] * v1d[1][1] * g1d[2][0];
    dq3d[0][2] = g1d[0][0] * v1d[1][2] * v1d[2][0];
    dq3d[1][2] = v1d[0][0] * g1d[1][2] * v1d[2][0];
    dq3d[2][2] = v1d[0][0] * v1d[1][2] * g1d[2][0];
    dq3d[0][3] = g1d[0][1] * v1d[1][0] * v1d[2][1];
    dq3d[1][3] = v1d[0][1] * g1d[1][0] * v1d[2][1];
    dq3d[2][3] = v1d[0][1] * v1d[1][0] * g1d[2][1];
    dq3d[0][4] = g1d[0][0] * v1d[1][1] * v1d[2][1];
    dq3d[1][4] = v1d[0][0] * g1d[1][1] * v1d[2][1];
    dq3d[2][4] = v1d[0][0] * v1d[1][1] * g1d[2][1];
    dq3d[0][5] = g1d[0][0] * v1d[1][0] * v1d[2][2];
    dq3d[1][5] = v1d[0][0] * g1d[1][0] * v1d[2][2];
    dq3d[2][5] = v1d[0][0] * v1d[1][0] * g1d[2][2];
}

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


__device__
void transform0(
    const int r, const int c, const int lj, const int li, 
    const double *cart, double *sphr) {
    switch (li) {
        case 0:
        case 1:
            switch (lj) {
                case 0:
                case 1:
                    // sphr = cart
                    memcpy(sphr, cart, r * c * sizeof(double));
                    break;
                case 2:
                    // sphr(1, :) = cart(3, :) - 0.5_wp * (cart(1, :) + cart(2, :))
                    for (int j = 0; j < c; j++) {
                        sphr[0 * c + j] = cart[2 * c + j] - 0.5 * (cart[0 * c + j] + cart[1 * c + j]);
                        sphr[1 * c + j] = S3 * cart[4 * c + j];
                        sphr[2 * c + j] = S3 * cart[5 * c + j];
                        sphr[3 * c + j] = S3_4 * (cart[0 * c + j] - cart[1 * c + j]);
                        sphr[4 * c + j] = S3 * cart[3 * c + j];
                    }
                    break;
                case 3:
                    // TODO: sphr = matmul(ftrafo, cart)
                    printf("[Error] Matmul with ftrafo not yet implemented.\n");
                    assert(false);
                    break;
                case 4:
                    // TODO: sphr = matmul(gtrafo, cart)
                    printf("[Error] Matmul with gtrafo not yet implemented.\n");
                    assert(false);
                    break;
                default:
                    printf("[Fatal] Moments higher than g are not supported.\n");
                    assert(false);
            }
            break;

        case 2:
            switch (lj) {
                case 0:
                case 1:
                    // sphr(:, 1) = cart(:, 3) - 0.5_wp * (cart(:, 1) + cart(:, 2))
                    // for (int i = 0; i < cols; i++) {
                    //     sphr[i * cols + 0] = cart[i * cols + 2] - 0.5 * (cart[i * cols + 0] + cart[i * cols + 1]);
                    //     sphr[i * cols + 1] = S3 * cart[i * cols + 4];
                    //     sphr[i * cols + 2] = S3 * cart[i * cols + 5];
                    //     sphr[i * cols + 3] = S3_4 * (cart[i * cols + 0] - cart[i * cols + 1]);
                    //     sphr[i * cols + 4] = S3 * cart[i * cols + 3];
                    // }
                    break;
                case 2:
                    // TODO: sphr = matmul(dtrafo, matmul(cart, transpose(dtrafo)))
                    printf("[Error] Matmul with dtrafo not yet implemented.\n");
                    assert(false);
                    break;
                case 3:
                    // TODO: sphr = matmul(ftrafo, matmul(cart, transpose(dtrafo)))
                    printf("[Error] Matmul with ftrafo and dtrafo not yet implemented.\n");
                    assert(false);
                    break;
                case 4:
                    // TODO: sphr = matmul(gtrafo, matmul(cart, transpose(dtrafo)))
                    printf("[Error] Matmul with gtrafo and dtrafo not yet implemented.\n");
                    assert(false);
                    break;
                default:
                    printf("[Fatal] Moments higher than g are not supported.\n");
                    assert(false);
            }
            break;

        case 3:
            switch (lj) {
                case 0:
                case 1:
                    // TODO: sphr = matmul(cart, transpose(ftrafo))
                    printf("[Error] Matmul with transpose(ftrafo) not yet implemented.\n");
                    assert(false);
                    break;
                case 2:
                    // TODO: sphr = matmul(dtrafo, matmul(cart, transpose(ftrafo)))
                    printf("[Error] Matmul with dtrafo and transpose(ftrafo) not yet implemented.\n");
                    assert(false);
                    break;
                case 3:
                    // TODO: sphr = matmul(ftrafo, matmul(cart, transpose(ftrafo)))
                    printf("[Error] Matmul with ftrafo and transpose(ftrafo) not yet implemented.\n");
                    assert(false);
                    break;
                case 4:
                    // TODO: sphr = matmul(gtrafo, matmul(cart, transpose(ftrafo)))
                    printf("[Error] Matmul with gtrafo and transpose(ftrafo) not yet implemented.\n");
                    assert(false);
                    break;
                default:
                    printf("[Fatal] Moments higher than g are not supported.\n");
                    assert(false);
            }
            break;

        case 4:
            switch (lj) {
                case 0:
                case 1:
                    // TODO: sphr = matmul(cart, transpose(gtrafo))
                    printf("[Error] Matmul with transpose(gtrafo) not yet implemented.\n");
                    assert(false);
                    break;
                case 2:
                    // TODO: sphr = matmul(dtrafo, matmul(cart, transpose(gtrafo)))
                    printf("[Error] Matmul with dtrafo and transpose(gtrafo) not yet implemented.\n");
                    assert(false);
                    break;
                case 3:
                    // TODO: sphr = matmul(ftrafo, matmul(cart, transpose(gtrafo)))
                    printf("[Error] Matmul with ftrafo and transpose(gtrafo) not yet implemented.\n");
                    assert(false);
                    break;
                case 4:
                    // TODO: sphr = matmul(gtrafo, matmul(cart, transpose(gtrafo)))
                    printf("[Error] Matmul with gtrafo and transpose(gtrafo) not yet implemented.\n");
                    assert(false);
                    break;
                default:
                    printf("[Fatal] Moments higher than g are not supported.\n");
                    assert(false);
            }
            break;

        default:
            printf("[Fatal] Moments higher than g are not supported.\n");
            assert(false);
    }
}

// pure subroutine transform1(lj, li, cart, sphr)
//    integer, intent(in) :: li
//    integer, intent(in) :: lj
//    real(wp), intent(in) :: cart(:, :, :)
//    real(wp), intent(out) :: sphr(:, :, :)
//    integer :: k

//    do k = 1, size(cart, 1)
//       call transform0(lj, li, cart(k, :, :), sphr(k, :, :))
//    end do
// end subroutine transform1

__device__
void transform1(
   const int r, const int c, const int d, const int lj, const int li, 
    const double *cart, double *sphr) {
    for (int k = 0; k < r; k++) {
        transform0(c, d, lj, li, &cart[k * c * d], &sphr[k * c * d]);
    }
}

// pure subroutine transform2(lj, li, cart, sphr)
//    integer, intent(in) :: li
//    integer, intent(in) :: lj
//    real(wp), intent(in) :: cart(:, :, :, :)
//    real(wp), intent(out) :: sphr(:, :, :, :)
//    integer :: k, l

//    do l = 1, size(cart, 2)
//       do k = 1, size(cart, 1)
//          call transform0(lj, li, cart(k, l, :, :), sphr(k, l, :, :))
//       end do
//    end do
// end subroutine transform2

__device__
void transform2(
    const int r, const int c, const int d, const int e, 
    const int lj, const int li, 
    const double *cart, double *sphr) {
    for (int l = 0; l < r; l++) {
        for (int k = 0; k < c; k++) {
            transform0(d, e, lj, li, &cart[(l * c + k) * d * e], &sphr[(l * c + k) * d * e]);
        }
    }
}