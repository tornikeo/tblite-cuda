#include "h0.h"
#include "integral/overlap.h"
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


// 2D array
// template <typename T, int A, int B>
// __device__ 
// void zero(T (&arr)[A][B])
// {
//   for (int i = 0; i < A; ++i)
//   {
//     for (int j = 0; j < B; ++j)
//     {
//       arr[i][j] = static_cast<T>(0);
//     }
//   }
// }

__device__ void get_hamiltonian(
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
  float stmp[msao(MAXL)][msao(MAXL)] = {0};
  float dtmpi[msao(MAXL)][msao(MAXL)][3] = {0}; // Notice the dims are reversed
  float qtmpi[msao(MAXL)][msao(MAXL)][6] = {0};

  zero(overlap);
  zero(dpint);
  zero(qpint);
  zero(hamiltonian);

  /*       do iat = 1, mol%nat
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
        call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
          & r2, vec, bas%intcut, stmp, dtmpi, qtmpi) */
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

      nao = msao(bas.cgto[jsh][jzp].ang);
      for (iao = 0; iao < msao(bas.cgto[ish][izp].ang); ++iao)
      {
        for (jao = 0; jao < nao; ++jao)
        {
          // ij = jao + nao * iao;
          
        }
      }
    }
  }
}