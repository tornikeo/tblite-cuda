#ifndef TBLITE_MAIN   /* Include guard */
#define TBLITE_MAIN
#define MAXL 6
#define MAXL2 12 // 2 x MAXL
#define SQRT_PI 1.77245385091  // Approximation of sqrt(pi)
#define SQRT_PI3 5.56832799683 // Approximation of sqrt(pi)^3
#define S3 sqrt(3.0)
#define S3_4 (S3 * 0.5)
#define D32 (3.0 / 2.0)
#define S3_8 sqrt(3.0 / 8.0)
#define S5_8 sqrt(5.0 / 8.0)
#define S6 sqrt(6.0)
#define S15 sqrt(15.0)
#define S15_4 sqrt(15.0 / 4.0)
#define S45 sqrt(45.0)
#define S45_8 sqrt(45.0 / 8.0)
#define D38 (3.0 / 8.0)
#define D34 (3.0 / 4.0)
#define S5_16 sqrt(5.0 / 16.0)
#define S10 sqrt(10.0)
#define S10_8 sqrt(10.0 / 8.0)
#define S35_4 sqrt(35.0 / 4.0)
#define S35_8 sqrt(35.0 / 8.0)
#define S35_64 sqrt(35.0 / 64.0)
#define S45_4 sqrt(45.0 / 4.0)
#define S315_8 sqrt(315.0 / 8.0)
#define S315_16 sqrt(315.0 / 16.0)
#pragma once

template<typename T, int Rows, int Cols>
class Matrix2D {
private:
    T data[Rows * Cols];

public:
    __device__ T& at(int row, int col) {
        return data[row * Cols + col];
    }

    __device__ const T& at(int row, int col) const {
        return data[row * Cols + col];
    }

    __device__ int rowCount() const { return Rows; }
    __device__ int colCount() const { return Cols; }

    // Optional raw access
    __device__ T* rawData() { return data; }
    __device__ const T* rawData() const { return data; }
};

template<typename T, int Depth, int Rows, int Cols>
class Matrix3D {
private:
    T data[Depth * Rows * Cols];

public:
    __device__ T& at(int depth, int row, int col) {
        return data[depth * Rows * Cols + row * Cols + col];
    }

    __device__ const T& at(int depth, int row, int col) const {
        return data[depth * Rows * Cols + row * Cols + col];
    }

    __device__ int depthCount() const { return Depth; }
    __device__ int rowCount() const { return Rows; }
    __device__ int colCount() const { return Cols; }

    // Optional raw access
    __device__ T* rawData() { return data; }
    __device__ const T* rawData() const { return data; }
};

template<typename T, int Time, int Depth, int Rows, int Cols>
class Matrix4D {
private:
    T data[Time * Depth * Rows * Cols];

public:
    __device__ T& at(int time, int depth, int row, int col) {
        return data[time * Depth * Rows * Cols + depth * Rows * Cols + row * Cols + col];
    }

    __device__ const T& at(int time, int depth, int row, int col) const {
        return data[time * Depth * Rows * Cols + depth * Rows * Cols + row * Cols + col];
    }

    __device__ int timeCount() const { return Time; }
    __device__ int depthCount() const { return Depth; }
    __device__ int rowCount() const { return Rows; }
    __device__ int colCount() const { return Cols; }

    // Optional raw access
    __device__ T* rawData() { return data; }
    __device__ const T* rawData() const { return data; }
};

/*
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
*/

#define MAXG 4
typedef struct {
    int ang = 0;
    int nprim = 0;
    double alpha[MAXG] = {0.0};
    double coeff[MAXG] = {0.0};
} cgto_type;


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

// #define MSAO 5
// #define MLAO 6
// #define LMAP 1

typedef struct {
    double* at;
    size_t R;
    size_t C;
} matrix;

typedef struct {
    double* at;
    size_t R;
    size_t C;
    size_t W;
} tens3;

typedef struct {
    double* at;
    size_t R;
    size_t C;
    size_t W;
    size_t H;
} tens4;
        
template <size_t N, size_t M, size_t MSAO, size_t MLAO>
__device__
void transform2(const int lj, const int li, 
    const double (&cart)[N][M][MLAO][MLAO], double (&sphr)[N][M][MSAO][MSAO]);
template <size_t N, size_t MSAO, size_t MLAO>
__device__
void transform1(const int lj, const int li, 
    const double (&cart)[N][MLAO][MLAO], double (&sphr)[N][MSAO][MSAO]);

// template <size_t msaoj, size_t msaoi, size_t mlaoi, size_t mlaoj>
__device__
void transform0(
    const int lj, const int li, 
    const matrix &cart, matrix &sphr);

// #define MSAO 3
// #define MLAO 3
// template <size_t MSAO>//, size_t MLAO, size_t LMAP>
// template <size_t MSAO>
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
//     double (&doverlap)[3][MSAO][MSAO],
//     double (&ddpinti)[3][3][MSAO][MSAO],
//     double (&dqpinti)[3][6][MSAO][MSAO],
//     double (&ddpintj)[3][3][MSAO][MSAO],
//     double (&dqpintj)[3][6][MSAO][MSAO]
// );

__global__ void test_call_multipole_grad_cgto();
#endif // FOO_H_