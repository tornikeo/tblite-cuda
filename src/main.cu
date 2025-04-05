// type :: adjacency_list
// !> Offset index in the neighbour map
// integer, allocatable :: inl(:)
// !> Number of neighbours for each atom
// integer, allocatable :: nnl(:)
// !> Index of the neighbouring atom
// integer, allocatable :: nlat(:)
// !> Cell index of the neighbouring atom
// integer, allocatable :: nltr(:)
// end type adjacency_list

typedef struct {
    int *inl;   // Offset index in the neighbour map
    int *nnl;   // Number of neighbours for each atom
    int *nlat;  // Index of the neighbouring atom
    int *nltr;  // Cell index of the neighbouring atom
} adjacency_list;


// !> Structure representation
// type :: structure_type

//    !> Number of atoms
//    integer :: nat = 0

//    !> Number of unique species
//    integer :: nid = 0

//    !> Number of bonds
//    integer :: nbd = 0

//    !> Species identifier
//    integer, allocatable :: id(:)

//    !> Atomic number for each species
//    integer, allocatable :: num(:)

//    !> Element symbol for each species
//    character(len=symbol_length), allocatable :: sym(:)

//    !> Cartesian coordinates, in Bohr
//    real(wp), allocatable :: xyz(:, :)

//    !> Number of unpaired electrons
//    integer :: uhf = 0

//    !> Total charge
//    real(wp) :: charge = 0.0_wp

//    !> Lattice parameters
//    real(wp), allocatable :: lattice(:, :)

//    !> Periodic directions
//    logical, allocatable :: periodic(:)

//    !> Bond indices
//    integer, allocatable :: bond(:, :)

// end type structure_type

typedef struct {
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


// integer, parameter :: maxl = 6
// integer, parameter :: maxl2 = maxl*2
// integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
// integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
// integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
// real(wp), parameter :: sqrtpi = sqrt(pi)
// real(wp), parameter :: sqrtpi3 = sqrtpi**3

#define MAXL 6
#define MAXL2 12 // 2 x MAXL
#define SQRT_PI 1.77245385091 // Approximation of sqrt(pi)
#define SQRT_PI3 5.56832799683 // Approximation of sqrt(pi)^3


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

__device__ void form_product(const double *a, const double *b, int la, int lb, double *d) {
    if (la >= 4 || lb >= 4) goto level_40;
    if (la >= 3 || lb >= 3) goto level_30;
    if (la >= 2 || lb >= 2) goto level_20;

    // <s|s> = <s>
    d[0] = a[0] * b[0];
    if (la == 0 && lb == 0) return;

    // <s|p> = <s|*(|s>+|p>) = <s> + <p>
    d[1] = a[0] * b[1] + a[1] * b[0];
    if (la == 0 || lb == 0) return;

    // <p|p> = (<s|+<p|)*(|s>+|p>) = <s> + <p> + <d>
    d[2] = a[1] * b[1];
    return;

level_20:
    // <s|d> = <s|*(|s>+|p>+|d>) = <s> + <p> + <d>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    if (la == 0 || lb == 0) return;

    // <p|d> = (<s|+<p|)*(|s>+|p>+|d>) = <s> + <p> + <d> + <f>
    d[2] += a[1] * b[1];
    d[3] = a[1] * b[2] + a[2] * b[1];
    if (la <= 1 || lb <= 1) return;

    // <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>) = <s> + <p> + <d> + <f> + <g>
    d[4] = a[2] * b[2];
    return;

level_30:
    // <s|f> = <s|*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f>
    d[0] = a[0] * b[0];
    d[1] = a[0] * b[1] + a[1] * b[0];
    d[2] = a[0] * b[2] + a[2] * b[0];
    d[3] = a[0] * b[3] + a[3] * b[0];
    if (la == 0 || lb == 0) return;

    // <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g>
    d[2] += a[1] * b[1];
    d[3] += a[1] * b[2] + a[2] * b[1];
    d[4] = a[1] * b[3] + a[3] * b[1];
    if (la <= 1 || lb <= 1) return;

    // <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g> + <h>
    d[4] += a[2] * b[2];
    d[5] = a[2] * b[3] + a[3] * b[2];
    if (la <= 2 || lb <= 2) return;

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
    if (la == 0 || lb == 0) return;

    // <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h>
    d[2] += a[1] * b[1];
    d[3] += a[1] * b[2] + a[2] * b[1];
    d[4] += a[1] * b[3] + a[3] * b[1];
    d[5] = a[1] * b[4] + a[4] * b[1];
    if (la <= 1 || lb <= 1) return;

    // <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
    d[4] += a[2] * b[2];
    d[5] += a[2] * b[3] + a[3] * b[2];
    d[6] = a[2] * b[4] + a[4] * b[2];
    if (la <= 2 || lb <= 2) return;

    // <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
    d[6] += a[3] * b[3];
    d[7] = a[3] * b[4] + a[4] * b[3];
    if (la <= 3 || lb <= 3) return;

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

__device__ void horizontal_shift(double ae, int l, double *cfs) {
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
        default:
            break;
    }
}

