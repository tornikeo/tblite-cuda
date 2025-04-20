#ifndef H_XTB
#define X_XTB

#define MAXG 12

typedef struct {
  int nat;
  int nid;
  int nbd;
  int *id;
  int *num;
  float *xyz;
} structure_type;

// integer, parameter :: maxg = 12
// !> Contracted Gaussian type basis function
// type :: cgto_type
//   !> Angular momentum of this basis function
//   integer :: ang = -1
//   !> Contraction length of this basis function
//   integer :: nprim = 0
//   !> Exponent of the primitive Gaussian functions
//   real(wp) :: alpha(maxg) = 0.0_wp
//   !> Contraction coefficients of the primitive Gaussian functions,
//   !> might contain normalization
//   real(wp) :: coeff(maxg) = 0.0_wp
// end type cgto_type
typedef struct {
  int ang;
  int nprim;
  float alpha[MAXG];
  float coeff[MAXG];
} cgto_type;

// type :: basis_type
//   !> Maximum angular momentum of all basis functions,
//   !> used to determine scratch size in integral calculation
//   integer :: maxl = 0
//   !> Number of shells in this basis set
//   integer :: nsh = 0
//   !> Number of spherical atomic orbitals in this basis set
//   integer :: nao = 0
//   !> Integral cutoff as maximum exponent of Gaussian product theoreom to consider
//   real(wp) :: intcut = 0.0_wp
//   !> Smallest primitive exponent in the basis set
//   real(wp) :: min_alpha = huge(0.0_wp)
//   !> Number of shells for each species
//   integer, allocatable :: nsh_id(:)
//   !> Number of shells for each atom
//   integer, allocatable :: nsh_at(:)
//   !> Number of spherical atomic orbitals for each shell
//   integer, allocatable :: nao_sh(:)
//   !> Index offset for each shell in the atomic orbital space
//   integer, allocatable :: iao_sh(:)
//   !> Index offset for each atom in the shell space
//   integer, allocatable :: ish_at(:)
//   !> Mapping from spherical atomic orbitals to the respective atom
//   integer, allocatable :: ao2at(:)
//   !> Mapping from spherical atomic orbitals to the respective shell
//   integer, allocatable :: ao2sh(:)
//   !> Mapping from shells to the respective atom
//   integer, allocatable :: sh2at(:)
//   !> Contracted Gaussian basis functions forming the basis set
//   type(cgto_type), allocatable :: cgto(:, :)
// end type basis_type
typedef struct 
{
  int maxl = 0;
  int nsh = 0;
  int nao = 0;
  float intcut = 0;
  // min_alpha don't care
  int *nsh_id;
  int *nsh_at; 
  int *nao_sh;
  int *iao_sh;
  int *ish_at;
  int *ao2at;
  int *ao2sh;
  int *sh2at;
  cgto_type *cgto;
} basis_type;

/*
   type :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq1(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq2(:, :)
      !> Enhancement factor to scale the Hamiltonian elements
      real(wp), allocatable :: hscale(:, :, :, :)
      !> Polynomial coefficients for distance dependent enhancement factor
      real(wp), allocatable :: shpoly(:, :)
      !> Atomic radius for polynomial enhancement
      real(wp), allocatable :: rad(:)
      !> Reference occupation numbers
      real(wp), allocatable :: refocc(:, :)
   end type tb_hamiltonian

   subroutine new_hamiltonian(self, mol, bas, spec)
      type(tb_hamiltonian), intent(out) :: self
      type(structure_type), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      class(tb_h0spec), intent(in) :: spec

      integer :: mshell

      mshell = maxval(bas%nsh_id)
      allocate(self%selfenergy(mshell, mol%nid), self%kcn(mshell, mol%nid), &
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
      call spec%get_selfenergy(mol, bas, self%selfenergy)
      call spec%get_cnshift(mol, bas, self%kcn)
      call spec%get_q1shift(mol, bas, self%kq1)
      call spec%get_q2shift(mol, bas, self%kq2)

      allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
      call spec%get_hscale(mol, bas, self%hscale)

      allocate(self%shpoly(mshell, mol%nid), self%rad(mol%nid))
      call spec%get_rad(mol, bas, self%rad)
      call spec%get_shpoly(mol, bas, self%shpoly)

      allocate(self%refocc(mshell, mol%nid))
      call spec%get_reference_occ(mol, bas, self%refocc)
   end subroutine new_hamiltonian
*/
typedef struct {
  float *selfenergy; // Flattened 2D array, [mshell, mol%nid]
  float *kcn;        // Flattened 2D array, [mshell, mol%nid]
  float *kq1;        // Flattened 2D array, [mshell, mol%nid]
  float *kq2;        // Flattened 2D array, [mshell, mol%nid]
  float *hscale;     // Flattened 4D array, [mshell, mshell, mol%nid, mol%nid]
  float *shpoly;     // Flattened 2D array, [mshell, mol%nid]
  float *rad;        // 1D array, [mol%nid]
  float *refocc;     // Flattened 2D array, [mshell, mol%nid]
} tb_hamiltonian;

typedef struct
{
  // type(basis_type) :: bas
  // type(tb_hamiltonian) :: h0
  // class(ncoord_type), allocatable :: ncoord
  // type(tb_repulsion), allocatable :: repulsion
  // type(tb_coulomb), allocatable :: coulomb
  // type(halogen_correction), allocatable :: halogen
  // class(dispersion_type), allocatable :: dispersion
  // real(wp) :: mixer_damping = mixer_damping_default
  // integer :: max_iter = max_iter_default
  // logical :: save_integrals = .false.
  // !> List of additional interaction containers
  // type(container_list), allocatable :: interactions
  basis_type bas;
  tb_hamiltonian h0;
} xtb_calculator;

typedef struct
{
  
} wavefunction_type;


#endif