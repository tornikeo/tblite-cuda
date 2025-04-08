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
    const int lj, const int li, 
    const matrix &cart, 
    matrix &sphr
) {
    switch (li) {
        case 0:
        case 1:
            switch (lj) {
                case 0:
                case 1:
                    // sphr = cart
                    // memcpy(sphr, cart, msaoj * msaoi * sizeof(double));
                    memcpy( sphr.at, cart.at, cart.R*cart.C*sizeof(double));
                    break;
                case 2:
                    // sphr(1, :) = cart(3, :) - 0.5_wp * (cart(1, :) + cart(2, :))
                    for (int i = 0; i < cart.C; i++) {
                        sphr.at[0*sphr.C+i] = cart.at[2*cart.C+i] - 0.5 * (cart.at[0*cart.C+i] + cart.at[1*cart.C+i]);
                        sphr.at[1*sphr.C+i] = S3 * cart.at[4*cart.C+i];
                        sphr.at[2*sphr.C+i] = S3 * cart.at[5*cart.C+i];
                        sphr.at[3*sphr.C+i] = S3_4 * (cart.at[0*cart.C+i] - cart.at[1*cart.C+i]);
                        sphr.at[4*sphr.C+i] = S3 * cart.at[3*cart.C+i];
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

template <size_t N, size_t MSAO, size_t MLAO>
__device__
void transform1(const int lj, const int li, 
    const double (&cart)[N][MLAO][MLAO], double (&sphr)[N][MSAO][MSAO]
) {
    // for (int i = 0; i < N; i++) {
    //     transform0(lj, li, cart[i], sphr[i]);
    // }
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

template <size_t N, size_t M, size_t MSAO, size_t MLAO>
__device__
void transform2(const int lj, const int li, 
    const double (&cart)[N][M][MLAO][MLAO], double (&sphr)[N][M][MSAO][MSAO]
) {
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < M; i++) {
    //         transform0(lj, li, cart[i][j], sphr[i][j]);
    //     }
    // }
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
double overlap_1d(const int moment, const double alpha) {
    const double dfactorial[] = {1.0, 1.0, 3.0, 15.0, 105.0, 945.0, 10395.0, 135135.0};
    if (moment % 2 == 0) {
        return pow(0.5 / alpha, moment / 2) * dfactorial[moment / 2];
    } else {
        return 0.0;
    }
}

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

constexpr size_t MSAO[7] = {1, 3, 5, 7, 9, 11, 13};
constexpr size_t MLAO[7] = {1, 3, 6, 10, 15, 21, 28};
constexpr size_t LMAP[7] = {0, 1, 4, 10, 20, 35, 56};
constexpr int LX[8][3] = {
    {0,0,0}, 
    {0,1,0},
    {0,0,1},
    {1,0,0},
    {2,0,0},
    {0,2,0},
    {0,0,2},
    {1,1,2},
};

template  <size_t msao_j, size_t msao_i>
__device__ 
void multipole_grad_cgto(
    const cgto_type cgtoj,
    const cgto_type cgtoi,
    const double r2, 
    const double vec[3],
    const double intcut,

    double (&overlap)[msao_j][msao_i],
    double (&dpint)[3][msao_j][msao_i],
    double (&qpint)[6][msao_j][msao_i],
    double (&doverlap)[3][msao_j][msao_i],
    double (&ddpinti)[3][3][msao_j][msao_i],
    double (&dqpinti)[3][6][msao_j][msao_i],
    double (&ddpintj)[3][3][msao_j][msao_i],
    double (&dqpintj)[3][6][msao_j][msao_i]
) {

    // double eab, oab, est, s1d[MAXL2 + 1], rpi[3], rpj[3], cc, val, dip[3], quad[6];
    // double pre, grad[3], ddip[3][3], dquad[3][6], tr, dtr[3];
    // double s3d[MLAO][MLAO] = {0.0};
    // double d3dj[3][MLAO][MLAO] = {0.0};
    // double q3dj[6][MLAO][MLAO] = {0.0};
    // double ds3d[3][MLAO][MLAO] = {0.0};
    // double dd3di[3][3][MLAO][MLAO] = {0.0};
    // double dq3di[3][6][MLAO][MLAO] = {0.0};
    // double dd3dj[3][3][MLAO][MLAO] = {0.0};
    // double dq3dj[3][6][MLAO][MLAO] = {0.0};

    // for (int ip = 0; ip < cgtoi.nprim; ip++) {
    //     for (int jp = 0; jp < cgtoj.nprim; jp++) {
    //         eab = cgtoi.alpha[ip] + cgtoj.alpha[jp];
    //         oab = 1.0 / eab;
    //         est = cgtoi.alpha[ip] * cgtoj.alpha[jp] * r2 * oab;
    //         if (est > intcut) continue;
    //         pre = exp(-est) * SQRT_PI3 * pow(sqrt(oab), 3);
    //         for (int i = 0; i < 3; i++) {
    //             rpi[i] = -vec[i] * cgtoj.alpha[jp] * oab;
    //             rpj[i] = vec[i] * cgtoi.alpha[ip] * oab;
    //         }
    //         for (int l = 0; l <= cgtoi.ang + cgtoj.ang + 3; l++) {
    //             s1d[l] = overlap_1d(l, eab);
    //         }
    //         cc = cgtoi.coeff[ip] * cgtoj.coeff[jp] * pre;
    //         for (int mli = 0; mli < MLAO; mli++) {
    //             for (int mlj = 0; mlj < MLAO; mlj++) {
    //                 multipole_grad_3d(
    //                     rpj, rpi, cgtoj.alpha[jp], cgtoi.alpha[ip],
    //                     lx[mlj + lmap[cgtoj.ang]], lx[mli + lmap[cgtoj.ang]],
    //                     s1d, val, dip, quad, grad, ddip, dquad
    //                 );
    //                 s3d[mlj][mli] += cc * val;
    //                 for (int i = 0; i < 3; i++) {
    //                     d3dj[i][mlj][mli] += cc * dip[i];
    //                     ds3d[i][mlj][mli] += cc * grad[i];
    //                     for (int j = 0; j < 3; j++) {
    //                         dd3dj[i][j][mlj][mli] += cc * ddip[i][j];
    //                     }
    //                     for (int j = 0; j < 6; j++) {
    //                         dq3dj[i][j][mlj][mli] += cc * dquad[i][j];
    //                     }
    //                 }
    //                 for (int i = 0; i < 6; i++) {
    //                     q3dj[i][mlj][mli] += cc * quad[i];
    //                 }
    //             }
    //         }
    //     }
    // }

    // // for (int mli = 0; mli < MLAO(cgtoi.ang); mli++) {
    // //     for (int mlj = 0; mlj < MLAO(cgtoj.ang); mlj++) {
    // //         shift_operator(
    // //             vec, s3d[mlj][mli], d3dj[:, mlj][mli], q3dj[:, mlj][mli],
    // //             ds3d[:, mlj][mli], dd3dj[:, :, mlj][mli], dq3dj[:, :, mlj][mli],
    // //             dd3di[:, :, mlj][mli], dq3di[:, :, mlj][mli]
    // //         );
    // //     }
    // // }

    // transform0(cgtoj.ang, cgtoi.ang, s3d, overlap);
    // transform1(cgtoj.ang, cgtoi.ang, d3dj, dpint);
    // transform1(cgtoj.ang, cgtoi.ang, q3dj, qpint);
    // transform1(cgtoj.ang, cgtoi.ang, ds3d, doverlap);
    // transform2(cgtoj.ang, cgtoi.ang, dd3dj, ddpintj);
    // transform2(cgtoj.ang, cgtoi.ang, dq3dj, dqpintj);
    // transform2(cgtoj.ang, cgtoi.ang, dd3di, ddpinti);
    // transform2(cgtoj.ang, cgtoi.ang, dq3di, dqpinti);

    // for (int mli = 0; mli < MSAO; mli++) {
    //     for (int mlj = 0; mlj < MSAO; mlj++) {
    //         tr = 0.5 * (qpint[0][mlj][mli] + qpint[2][mlj][mli] + qpint[5][mlj][mli]);
    //         qpint[0][mlj][mli] = 1.5 * qpint[0][mlj][mli] - tr;
    //         qpint[1][mlj][mli] = 1.5 * qpint[1][mlj][mli];
    //         qpint[2][mlj][mli] = 1.5 * qpint[2][mlj][mli] - tr;
    //         qpint[3][mlj][mli] = 1.5 * qpint[3][mlj][mli];
    //         qpint[4][mlj][mli] = 1.5 * qpint[4][mlj][mli];
    //         qpint[5][mlj][mli] = 1.5 * qpint[5][mlj][mli] - tr;

    //         for (int i = 0; i < 3; i++) {
    //             dtr[i] = dqpinti[i][0][mlj][mli] + dqpinti[i][2][mlj][mli] + dqpinti[i][5][mlj][mli];
    //             dqpinti[i][0][mlj][mli] = 1.5 * dqpinti[i][0][mlj][mli] - 0.5 * dtr[i];
    //             dqpinti[i][1][mlj][mli] = 1.5 * dqpinti[i][1][mlj][mli];
    //             dqpinti[i][2][mlj][mli] = 1.5 * dqpinti[i][2][mlj][mli] - 0.5 * dtr[i];
    //             dqpinti[i][3][mlj][mli] = 1.5 * dqpinti[i][3][mlj][mli];
    //             dqpinti[i][4][mlj][mli] = 1.5 * dqpinti[i][4][mlj][mli];
    //             dqpinti[i][5][mlj][mli] = 1.5 * dqpinti[i][5][mlj][mli] - 0.5 * dtr[i];

    //             dtr[i] = dqpintj[i][0][mlj][mli] + dqpintj[i][2][mlj][mli] + dqpintj[i][5][mlj][mli];
    //             dqpintj[i][0][mlj][mli] = 1.5 * dqpintj[i][0][mlj][mli] - 0.5 * dtr[i];
    //             dqpintj[i][1][mlj][mli] = 1.5 * dqpintj[i][1][mlj][mli];
    //             dqpintj[i][2][mlj][mli] = 1.5 * dqpintj[i][2][mlj][mli] - 0.5 * dtr[i];
    //             dqpintj[i][3][mlj][mli] = 1.5 * dqpintj[i][3][mlj][mli];
    //             dqpintj[i][4][mlj][mli] = 1.5 * dqpintj[i][4][mlj][mli];
    //             dqpintj[i][5][mlj][mli] = 1.5 * dqpintj[i][5][mlj][mli] - 0.5 * dtr[i];
    //         }
    //     }
    // }
}


__global__
void test_call_multipole_grad_cgto(void) {
    // Example test call for multipole_grad_cgto
    const int MSAO = 3;
    // const int MLAO = 3;
    const int LMAP = 1;

    cgto_type cgtoj = {/* Initialize with appropriate values */};
    cgto_type cgtoi = {/* Initialize with appropriate values */};
    double r2 = 1.0;
    double vec[3] = {0.1, 0.2, 0.3};
    double intcut = 1e-6;

    double overlap[MSAO][MSAO] = {0.0};
    double dpint[3][MSAO][MSAO] = {0.0};
    double qpint[6][MSAO][MSAO] = {0.0};
    double doverlap[3][MSAO][MSAO] = {0.0};
    double ddpinti[3][3][MSAO][MSAO] = {0.0};
    double dqpinti[3][6][MSAO][MSAO] = {0.0};
    double ddpintj[3][3][MSAO][MSAO] = {0.0};
    double dqpintj[3][6][MSAO][MSAO] = {0.0};

    multipole_grad_cgto(
        cgtoj, cgtoi, 
        r2, vec, intcut,
        overlap, dpint, qpint, doverlap, ddpinti, dqpinti, ddpintj, dqpintj
    );

    // Print or validate results as needed
    printf("Test call completed.\n");
}
// __device__ 
// void multipole_grad_cgto(
//     const cgto_type cgtoj,
//     const cgto_type cgtoi,
//     const double r2, 
//     const double vec[3],
//     const double intcut,

//     double (&overlap)[MSAO][MSAO],
//     double (&dpint)[3][MSAO][MSAO],
//     double (&qpint)[6][MSAO][MSAO],
// ) {
    
// }