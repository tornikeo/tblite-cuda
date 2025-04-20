module kernel
  use tblite_basis_type, only: wp
  implicit none
contains

attributes(global) subroutine test(a)
  integer, device :: a(*)
  integer :: i
  i = threadIdx%x
  a(i) = i
  return
end subroutine test

  
end module kernel