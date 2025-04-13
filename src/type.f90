module tblite_basis_type

  implicit none
     !> Single precision real numbers
  integer, parameter :: sp = selected_real_kind(6)

  !> Double precision real numbers
  integer, parameter :: dp = selected_real_kind(15)

  !> Wanted precision
  integer, parameter :: wp = dp

  integer, parameter :: maxg = 12

  type :: cgto_type
    !> Angular momentum of this basis function
    integer :: ang = -1
    !> Contraction length of this basis function
    integer :: nprim = 0
    !> Exponent of the primitive Gaussian functions
    real(wp) :: alpha(maxg) = 0.0_wp
    !> Contraction coefficients of the primitive Gaussian functions,
    !> might contain normalization
    real(wp) :: coeff(maxg) = 0.0_wp
  end type cgto_type

  !> Collection of information regarding the basis set of a system
  type :: basis_type
    !> Maximum angular momentum of all basis functions,
    !> used to determine scratch size in integral calculation
    integer :: maxl = 0
    !> Number of shells in this basis set
    integer :: nsh = 0
    !> Number of spherical atomic orbitals in this basis set
    integer :: nao = 0
    !> Integral cutoff as maximum exponent of Gaussian product theoreom to consider
    real(wp) :: intcut = 0.0_wp
    !> Smallest primitive exponent in the basis set
    real(wp) :: min_alpha = huge(0.0_wp)
    !> Number of shells for each species
    integer, allocatable :: nsh_id(:)
    !> Number of shells for each atom
    integer, allocatable :: nsh_at(:)
    !> Number of spherical atomic orbitals for each shell
    integer, allocatable :: nao_sh(:)
    !> Index offset for each shell in the atomic orbital space
    integer, allocatable :: iao_sh(:)
    !> Index offset for each atom in the shell space
    integer, allocatable :: ish_at(:)
    !> Mapping from spherical atomic orbitals to the respective atom
    integer, allocatable :: ao2at(:)
    !> Mapping from spherical atomic orbitals to the respective shell
    integer, allocatable :: ao2sh(:)
    !> Mapping from shells to the respective atom
    integer, allocatable :: sh2at(:)
    !> Contracted Gaussian basis functions forming the basis set
    type(cgto_type), allocatable :: cgto(:, :)
  end type basis_type

contains

end module tblite_basis_type