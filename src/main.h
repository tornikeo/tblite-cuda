#ifndef TBLITE_MAIN   /* Include guard */
#define TBLITE_MAIN
#define MAXL 6
#define MAXL2 12 // 2 x MAXL
#define SQRT_PI 1.77245385091  // Approximation of sqrt(pi)
#define SQRT_PI3 5.56832799683 // Approximation of sqrt(pi)^3

/* type :: adjacency_list
 !> Offset index in the neighbour map
 integer, allocatable :: inl(:)
 !> Number of neighbours for each atom
 integer, allocatable :: nnl(:)
 !> Index of the neighbouring atom
 integer, allocatable :: nlat(:)
 !> Cell index of the neighbouring atom
 integer, allocatable :: nltr(:)
 end type adjacency_list */

typedef struct
{
    int *inl;  // Offset index in the neighbour map
    int *nnl;  // Number of neighbours for each atom
    int *nlat; // Index of the neighbouring atom
    int *nltr; // Cell index of the neighbouring atom
} adjacency_list;


/* !> Structure representation
 type :: structure_type
    !> Number of atoms
    integer :: nat = 0
    !> Number of unique species
    integer :: nid = 0
    !> Number of bonds
    integer :: nbd = 0
    !> Species identifier
    integer, allocatable :: id(:)
    !> Atomic number for each species
    integer, allocatable :: num(:)
    !> Element symbol for each species
    character(len=symbol_length), allocatable :: sym(:)
    !> Cartesian coordinates, in Bohr
    real(wp), allocatable :: xyz(:, :)
    !> Number of unpaired electrons
    integer :: uhf = 0
    !> Total charge
    real(wp) :: charge = 0.0_wp
    !> Lattice parameters
    real(wp), allocatable :: lattice(:, :)
    !> Periodic directions
    logical, allocatable :: periodic(:)
    !> Bond indices
    integer, allocatable :: bond(:, :)
 end type structure_type */

typedef struct
{
    int nat;          // Number of atoms
    int nid;          // Number of unique species
    int nbd;          // Number of bonds
    int *id;          // Species identifier
    int *num;         // Atomic number for each species
    char **sym;       // Element symbol for each species
    double **xyz;     // Cartesian coordinates, in Bohr
    int uhf;          // Number of unpaired electrons
    double charge;    // Total charge
    double **lattice; // Lattice parameters
    int *periodic;    // Periodic directions (logical)
    int **bond;       // Bond indices
} structure_type;

/* integer, parameter :: maxl = 6
 integer, parameter :: maxl2 = maxl*2
 integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
 integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
 integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
 real(wp), parameter :: sqrtpi = sqrt(pi)
 real(wp), parameter :: sqrtpi3 = sqrtpi**3*/

__device__ void form_product(const double *a, const double *b, 
    int la, int lb, double *d);
__device__ void horizontal_shift(double ae, int l, double *cfs);
__device__ void shift_operator(const double *vec, double s, 
    const double *di, const double *qi, const double *ds, 
    const double ddi[3][3], const double dqi[3][6], 
    double ddj[3][3], double dqj[3][6]);
__device__ void multipole_grad_3d(
        const double rpi[3], const double rpj[3],
        const double ai, const double aj, 
        const int li[3], const int lj[3], const double s1d[MAXL2],
        
        double &s3d, double d3d[3], double q3d[3], double ds3d[3], 
        double dd3d[3][3], double dq3d[3][6]);

#endif // FOO_H_