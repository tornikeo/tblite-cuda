#include "h0.h"
#include "integral/overlap.h"
#include "integral/multipole.h"
#include <stdio.h>
// #include "utils/array.h" // TODO: NOT WORKING

/* UTILS */

// 2D array
template <typename T, size_t A, size_t B>
__device__ 
inline void zero(T (&arr)[A][B])
{
  for (size_t i = 0; i < A; ++i)
  {
    for (size_t j = 0; j < B; ++j)
    {
      arr[i][j] = static_cast<T>(0);
    }
  }
}

// 3D array
template <typename T, size_t A, size_t B, size_t C>
__device__ 
inline void zero(T (&arr)[A][B][C])
{
  for (size_t i = 0; i < A; ++i)
  {
    for (size_t j = 0; j < B; ++j)
    {
      for (size_t k = 0; k < C; ++k)
      {
        arr[i][j][k] = static_cast<T>(0);
      }
    }
  }
}

// 4D array (already implemented)
template <typename T, size_t A, size_t B, size_t C, size_t D>
__device__ 
inline void zero(T (&arr)[A][B][C][D])
{
  for (size_t i = 0; i < A; ++i)
  {
    for (size_t j = 0; j < B; ++j)
    {
      for (size_t k = 0; k < C; ++k)
      {
        for (size_t l = 0; l < D; ++l)
        {
          arr[i][j][k][l] = static_cast<T>(0);
        }
      }
    }
  }
}


/*subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Occupation number
      real(wp), intent(out) :: nocc
      !> Reference occupation for each atom
      real(wp), intent(out) :: n0at(:)
      !> Reference occupation for each shell
      real(wp), intent(out) :: n0sh(:)

      integer :: iat, ish, izp, ii

      nocc = -mol%charge
      n0at(:) = 0.0_wp
      n0sh(:) = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            nocc = nocc + h0%refocc(ish, izp)
            n0at(iat) = n0at(iat) + h0%refocc(ish, izp)
            n0sh(ii+ish) = n0sh(ii+ish) + h0%refocc(ish, izp)
         end do
      end do

   end subroutine get_occupation
   */
__device__ void get_occupation(
    const structure_type &mol,
    const basis_type &bas,
    const tb_hamiltonian &h0,
    float &nocc,
    float (&n0at)[MAX_NAT],
    float (&n0sh)[MAX_NSH])
{
  int iat, ish, izp, ii;

  nocc = -mol.charge;
  for (int i = 0; i < MAX_NAT; ++i)
  {
    n0at[i] = 0.0f;
  }
  for (int i = 0; i < MAX_NSH; ++i)
  {
    n0sh[i] = 0.0f;
  }

  for (iat = 0; iat < mol.nat; ++iat)
  {
    izp = mol.id[iat];
    ii = bas.ish_at[iat];
    for (ish = 0; ish < bas.nsh_id[izp]; ++ish)
    {
      nocc += h0.refocc[ish][izp];
      n0at[iat] += h0.refocc[ish][izp];
      n0sh[ii + ish] += h0.refocc[ish][izp];
    }
  }
}

__device__ void get_selfenergy(
    const tb_hamiltonian &h0,
    const int (&id)[MAX_NAT],
    const int (&ish_at)[MAX_NAT],
    const int (&nshell)[MAX_NELEM],
    const float (&cn)[MAX_NAT],
    // const float (&qat)[MAX_NAT],
    float selfenergy[MAX_NSH], // static size
    float dsedcn[MAX_NSH]      // static size
                               // float dsedq[MAX_NSH],      // static size
)
{
  /* subroutine get_selfenergy(h0, id, ish_at, nshell, cn, qat, selfenergy, dsedcn, dsedq)
      type(tb_hamiltonian), intent(in) :: h0
      integer, intent(in) :: id(:)
      integer, intent(in) :: ish_at(:)
      integer, intent(in) :: nshell(:)
      real(wp), intent(in), optional :: cn(:)
      real(wp), intent(in), optional :: qat(:)
      real(wp), intent(out) :: selfenergy(:)
      real(wp), intent(out), optional :: dsedcn(:)
      real(wp), intent(out), optional :: dsedq(:)

      integer :: iat, izp, ish, ii

      selfenergy(:) = 0.0_wp
      if (present(dsedcn)) dsedcn(:) = 0.0_wp
      if (present(dsedq)) dsedq(:) = 0.0_wp
      do iat = 1, size(id)
         izp = id(iat)
         ii = ish_at(iat)
         do ish = 1, nshell(izp)
            selfenergy(ii+ish) = h0%selfenergy(ish, izp)
         end do
      end do
      if (present(cn)) then
         if (present(dsedcn)) then
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
                  dsedcn(ii+ish) = -h0%kcn(ish, izp)
               end do
            end do
         else */
  for (int i = 0; i < MAX_NSH; ++i)
  {
    selfenergy[i] = 0.0f;
    dsedcn[i] = 0.0f;
  }

  for (int iat = 0; iat < MAX_NAT; ++iat)
  {
    int izp = id[iat];
    int ii = ish_at[iat];
    for (int ish = 0; ish < nshell[izp]; ++ish)
    {
      selfenergy[ii + ish] = h0.selfenergy[ish][izp];
    }
  }

  for (int iat = 0; iat < MAX_NAT; ++iat)
  {
    int izp = id[iat];
    int ii = ish_at[iat];
    for (int ish = 0; ish < nshell[izp]; ++ish)
    {
      selfenergy[ii + ish] -= h0.kcn[ish][izp] * cn[iat];
      dsedcn[ii + ish] = -h0.kcn[ish][izp];
    }
  }
  /* TODO: Low priority */
  /* Other options qat, dsedq are not implemented */
}

/* 
!> Shift multipole operator from Ket function (center i) to Bra function (center j),
!> the multipole operator on the Bra function can be assembled from the lower moments
!> on the Ket function and the displacement vector using horizontal shift rules.
   pure subroutine shift_operator(vec, s, di, qi, dj, qj)
      !> Displacement vector of center i and j
      real(wp),intent(in) :: vec(:)
      !> Overlap integral between basis functions
      real(wp),intent(in) :: s
      !> Dipole integral with operator on Ket function (center i)
      real(wp),intent(in) :: di(:)
      !> Quadrupole integral with operator on Ket function (center i)
      real(wp),intent(in) :: qi(:)
      !> Dipole integral with operator on Bra function (center j)
      real(wp),intent(out) :: dj(:)
      !> Quadrupole integral with operator on Bra function (center j)
      real(wp),intent(out) :: qj(:)

      real(wp) :: tr

      ! Create dipole operator on Bra function from Ket function and shift contribution
      ! due to monopol displacement
      dj(1) = di(1) + vec(1)*s
      dj(2) = di(2) + vec(2)*s
      dj(3) = di(3) + vec(3)*s

      ! For the quadrupole operator on the Bra function we first construct the shift
      ! contribution from the dipole and monopol displacement, since we have to remove
      ! the trace contribution from the shift and the moment integral on the Ket function
      ! is already traceless
      qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
      qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
      qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
      qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
      qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
      qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
      ! Now collect the trace of the shift contribution
      tr = 0.5_wp * (qj(1) + qj(3) + qj(6))

      ! Finally, assemble the quadrupole operator on the Bra function from the operator
      ! on the Ket function and the traceless shift contribution
      qj(1) = qi(1) + 1.5_wp * qj(1) - tr
      qj(2) = qi(2) + 1.5_wp * qj(2)
      qj(3) = qi(3) + 1.5_wp * qj(3) - tr
      qj(4) = qi(4) + 1.5_wp * qj(4)
      qj(5) = qi(5) + 1.5_wp * qj(5)
      qj(6) = qi(6) + 1.5_wp * qj(6) - tr
   end subroutine shift_operator

*/
__device__
inline void shift_operator(
  const float (&vec)[3],
  const float s,
  const float (&di)[3],
  const float (&qi)[6],
  float (&dj)[3],
  float (&qj)[6]
)
{
  float tr;

  // Create dipole operator on Bra function from Ket function and shift contribution
  // due to monopole displacement
  dj[0] = di[0] + vec[0] * s;
  dj[1] = di[1] + vec[1] * s;
  dj[2] = di[2] + vec[2] * s;

  // For the quadrupole operator on the Bra function, construct the shift contribution
  // from the dipole and monopole displacement
  qj[0] = 2.0f * vec[0] * di[0] + vec[0] * vec[0] * s;
  qj[2] = 2.0f * vec[1] * di[1] + vec[1] * vec[1] * s;
  qj[5] = 2.0f * vec[2] * di[2] + vec[2] * vec[2] * s;
  qj[1] = vec[0] * di[1] + vec[1] * di[0] + vec[0] * vec[1] * s;
  qj[3] = vec[0] * di[2] + vec[2] * di[0] + vec[0] * vec[2] * s;
  qj[4] = vec[1] * di[2] + vec[2] * di[1] + vec[1] * vec[2] * s;

  // Collect the trace of the shift contribution
  tr = 0.5f * (qj[0] + qj[2] + qj[5]);

  // Assemble the quadrupole operator on the Bra function from the operator
  // on the Ket function and the traceless shift contribution
  qj[0] = qi[0] + 1.5f * qj[0] - tr;
  qj[1] = qi[1] + 1.5f * qj[1];
  qj[2] = qi[2] + 1.5f * qj[2] - tr;
  qj[3] = qi[3] + 1.5f * qj[3];
  qj[4] = qi[4] + 1.5f * qj[4];
  qj[5] = qi[5] + 1.5f * qj[5] - tr;
}

// NOTE: This function is an abomination and get_hamiltonian_gradient is even worse.
__device__ 
void get_hamiltonian(
    const structure_type &mol,
    /* const float (&trans)[MAX_TRANS][3],*/ // Unimplemented
    const adjacency_list &alist,
    const basis_type &bas,
    const tb_hamiltonian &h0,
    const float (&selfenergy)[MAX_NSH],
    float (&overlap)[MAX_NAO][MAX_NAO],
    float (&dpint)[MAX_NAO][MAX_NAO][3],
    float (&qpint)[MAX_NAO][MAX_NAO][6],
    float (&hamiltonian)[MAX_NAO][MAX_NAO])
{
  // printf("this is going to hurt.\n");
  /*
    integer :: i,j,l;
    integer :: iat, jat, izp, jzp, itr, k, img, inl
    integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
    real(wp) :: rr, r2, vec(3), cutoff2, hij, shpoly, dtmpj(3), qtmpj(6)
    real(wp), allocatable :: stmp(:), dtmpi(:, :), qtmpi(:, :)

    overlap(:, :) = 0.0_wp
    dpint(:, :, :) = 0.0_wp
    qpint(:, :, :) = 0.0_wp
    hamiltonian(:, :) = 0.0_wp

    */
  int i, j, l;
  int iat, jat, izp, jzp, itr, k, img, inl;
  int ish, jsh, is, js, ii, jj, iao, jao, nao, ij;
  float rr, r2, vec[3], cutoff2, hij, shpoly, dtmpj[3], qtmpj[6];
   
  /*
  allocate(stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2), qtmpi(6, msao(bas%maxl)**2))
  */ 
  float stmp[msao[MAXL]][msao[MAXL]] = {0};
  float dtmpi[msao[MAXL]][msao[MAXL]][3] = {0}; // Notice the dims are reversed
  float qtmpi[msao[MAXL]][msao[MAXL]][6] = {0};

  zero(overlap);
  zero(dpint);
  zero(qpint);
  zero(hamiltonian);

  /*       
  do iat = 1, mol%nat
  izp = mol%id(iat)
  is = bas%ish_at(iat)
  inl = alist%inl(iat)
  do img = 1, alist%nnl(iat)
    jat = alist%nlat(img+inl)
    itr = alist%nltr(img+inl)
    jzp = mol%id(jat)
    js = bas%ish_at(jat)
  */
  for (iat = 0; iat < mol.nat; ++iat)
  {
    izp = mol.id[iat];
    is = bas.ish_at[iat];
    inl = alist.inl[iat];
    for (img = 0; img < alist.nnl[iat]; ++img)
    {
      jat = alist.nlat[img + inl];
      itr = alist.nltr[img + inl];
      jzp = mol.id[jat];
      js = bas.ish_at[jat];
      /*
      vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
      ! write (*,*) "FORTRAN: our own vec", vec

      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
      */
      for (int d = 0; d < 3; ++d)
      {
        vec[d] = mol.xyz[d][iat] - mol.xyz[d][jat]; //TODO - alist.trans[d][itr]; lattice unimplemented
      }

      r2 = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
      rr = sqrt(sqrt(r2) / (h0.rad[jzp] + h0.rad[izp]));
      
      /* 
      do ish = 1, bas%nsh_id(izp)
        ii = bas%iao_sh(is+ish)
        do jsh = 1, bas%nsh_id(jzp)
          jj = bas%iao_sh(js+jsh)
      */
      for (int ish = 0; ish < bas.nsh_id[izp]; ish++)
      {
        ii = bas.iao_sh[is + ish];
        for (int jsh = 0; jsh < bas.nsh_id[jzp]; jsh++)
        {
          jj = bas.iao_sh[js + jsh];
          /* 
            call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
              & r2, vec, bas%intcut, stmp, dtmpi, qtmpi) */
          multipole_cgto(bas.cgto[jsh][jzp], bas.cgto[ish][izp], 
            r2, vec, bas.intcut, stmp, dtmpi, qtmpi);
          /* 
            shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                * (1.0_wp + h0%shpoly(jsh, jzp)*rr)

            hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                * h0%hscale(jsh, ish, jzp, izp) * shpoly

            nao = msao(bas%cgto(jsh, jzp)%ang)
            do iao = 1, msao(bas%cgto(ish, izp)%ang)
                do jao = 1, nao
                  ij = jao + nao*(iao-1)
          */
          shpoly = (1.0f + h0.shpoly[ish][izp] * rr) *
              (1.0f + h0.shpoly[jsh][jzp] * rr);

          hij = 0.5f * (selfenergy[is + ish] + selfenergy[js + jsh]) *
            h0.hscale[jsh][ish][jzp][izp] * shpoly;

          nao = msao[bas.cgto[jsh][jzp].ang];
          for (iao = 0; iao < msao[bas.cgto[ish][izp].ang]; ++iao)
          {
            for (jao = 0; jao < nao; ++jao)
            {
              // ij = jao + nao * iao;
              /* TODO later
                call shift_operator(vec, stmp(ij), dtmpi(:, ij), qtmpi(:, ij), &
                            & dtmpj, qtmpj)
              */
              shift_operator(vec, stmp[iao][jao], dtmpi[iao][jao], qtmpi[iao][jao],
                dtmpj, qtmpj
              );
              /*
              overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                  + stmp(ij)

              do k = 1, 3
                  ! $omp atomic
                  dpint(k, jj+jao, ii+iao) = dpint(k, jj+jao, ii+iao) &
                    + dtmpi(k, ij)
              end do

              do k = 1, 6
                  ! $omp atomic
                  qpint(k, jj+jao, ii+iao) = qpint(k, jj+jao, ii+iao) &
                    + qtmpi(k, ij)
              end do

              ! $omp atomic
              hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                  + stmp(ij) * hij

              */

              overlap[jj + jao][ii + iao] += stmp[iao][jao];

              for (k = 0; k < 3; ++k)
              {
                dpint[ii + iao][jj + jao][k] += dtmpi[iao][jao][k];
              }

              for (k = 0; k < 6; ++k)
              {
                qpint[ii + iao][jj + jao][k] += qtmpi[iao][jao][k];
              }

              hamiltonian[ii + iao][jj + jao] += stmp[iao][jao] * hij;

              /* 
              if (iat /= jat) then
                  ! $omp atomic
                  overlap(ii+iao, jj+jao) = overlap(ii+iao, jj+jao) &
                    + stmp(ij)
                  do k = 1, 3
                    ! $omp atomic
                    dpint(k, ii+iao, jj+jao) = dpint(k, ii+iao, jj+jao) &
                        + dtmpj(k)
                  end do

                  do k = 1, 6
                    ! $omp atomic
                    qpint(k, ii+iao, jj+jao) = qpint(k, ii+iao, jj+jao) &
                        + qtmpj(k)
                  end do
                  ! $omp atomic
                  hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                    + stmp(ij) * hij
              end if
              */
              if (iat != jat)
              {
                overlap[ii + iao][jj + jao] += stmp[iao][jao];

                for (k = 0; k < 3; ++k)
                {
                  dpint[ii + iao][jj + jao][k] += dtmpj[k];
                }

                for (k = 0; k < 6; ++k)
                {
                  qpint[ii + iao][jj + jao][k] += qtmpj[k];
                }

                hamiltonian[ii + iao][jj + jao] += stmp[iao][jao] * hij;
              }
            }
          }

              
        }
      }
    }
  }
}

__global__
void test_hamiltonian()
{

}