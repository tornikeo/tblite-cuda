program main
  use tblite_integral_multipole, only: multipole_grad_cgto
  use tblite_basis_type, only: cgto_type, wp
  use kernel, only: test
  use cudafor
  implicit none

  type(cgto_type), device :: cgtoi, cgtoj
  real(wp) :: r2, intcut
  real(wp), dimension(3), device :: vec
  real(wp), allocatable, device :: overlap(:,:), dpint(:,:,:), qpint(:,:,:)
  real(wp), allocatable, device :: doverlap(:,:,:), ddpinti(:,:,:,:), dqpinti(:,:,:,:)
  real(wp), allocatable, device :: ddpintj(:,:,:,:), dqpintj(:,:,:,:)
  integer :: msao_i, msao_j

  ! integer :: i
  ! integer, parameter :: n = 100
  ! integer, allocatable, device :: iarr(:)
  ! integer, allocatable :: iarr_h(:)

  ! allocate(iarr(n))
  ! allocate(iarr_h(n))

  ! iarr_h(:) = 0
  ! call test<<<1,n>>>(iarr)
  ! iarr_h = iarr

  ! print*, &
  ! "Errors: ", count(iarr_h.ne.(/ (i, i=1,n) /))
  
  ! print*, iarr_h
  ! deallocate(iarr)

  ! Initialize dummy data
  cgtoi%ang = 1
  cgtoi%nprim = 1
  cgtoi%alpha(1) = 1.0_wp
  cgtoi%coeff(1) = 1.0_wp

  cgtoj%ang = 1
  cgtoj%nprim = 1
  cgtoj%alpha(1) = 1.0_wp
  cgtoj%coeff(1) = 1.0_wp

  r2 = 1.0_wp
  vec = [1.0_wp, 0.0_wp, 0.0_wp]
  intcut = 10.0_wp

  msao_i = 3  ! Number of spherical atomic orbitals for cgtoi
  msao_j = 3  ! Number of spherical atomic orbitals for cgtoj

  allocate(overlap(msao_j, msao_i))
  allocate(dpint(3, msao_j, msao_i))
  allocate(qpint(6, msao_j, msao_i))
  allocate(doverlap(3, msao_j, msao_i))
  allocate(ddpinti(3, 3, msao_j, msao_i))
  allocate(dqpinti(3, 6, msao_j, msao_i))
  allocate(ddpintj(3, 3, msao_j, msao_i))
  allocate(dqpintj(3, 6, msao_j, msao_i))

  ! Call the subroutine
  call <<< 1,1 >>> multipole_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint, &
       & doverlap, ddpintj, dqpintj, ddpinti, dqpinti)

  ! Print some results
  print*, "Overlap integrals:"
  print*, overlap

  print*, "Dipole moment integrals:"
  print*, dpint

  print*, "Quadrupole moment integrals:"
  print*, qpint

  ! Deallocate arrays
  deallocate(overlap, dpint, qpint, doverlap, ddpinti, dqpinti, ddpintj, dqpintj)

end program main