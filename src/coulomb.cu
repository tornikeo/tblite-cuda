
#include "coulomb.h"

__device__
void get_amat_0d(
  const structure_type &mol,
  const int (&nshell)[MAX_NAT], 
  const int (&offset)[MAX_NSH], 
  const float (&hubbard)[MSHELL][MSHELL][MAX_NELEM][MAX_NELEM], 
  const float gexp,
  float (&amat)[MAX_NSH][MAX_NSH]
) {
  /*do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         do ish = 1, nshell(iat)
            do jsh = 1, nshell(jat)
               gam = hubbard(jsh, ish, jzp, izp)
               tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
               ! $omp atomic
               amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + tmp
               ! $omp atomic
               amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + tmp
            end do
         end do
      end do
      do ish = 1, nshell(iat)
         do jsh = 1, ish-1
            gam = hubbard(jsh, ish, izp, izp)
            ! $omp atomic
            amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + gam
            ! $omp atomic
            amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + gam
         end do
         ! $omp atomic
         amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + hubbard(ish, ish, izp, izp)
      end do
   end do
  */
  int iat, jat, izp, jzp, ii, jj, ish, jsh;
  float vec[3], r1, r1g, gam, tmp;

  for (int iat = 0; iat < mol.nat; iat++)
  {
    for (int jat = 0; jat < iat; jat++) {
      int jzp = mol.id[jat];
      int jj = offset[jat];
      for (int k = 0; k < 3; k++) {
        vec[k] = mol.xyz[jat][k] - mol.xyz[iat][k];
      }
      r1 = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      r1g = powf(r1, gexp);
      for (int ish = 0; ish < nshell[iat]; ish++) {
        for (int jsh = 0; jsh < nshell[jat]; jsh++) {
          gam = hubbard[ish][jsh][izp][jzp];
          tmp = 1.0f / powf(r1g + powf(gam, -gexp), 1.0f / gexp);
          amat[ii + ish][jj + jsh] += tmp;
          amat[jj + jsh][ii + ish] += tmp;
        }
      }
    }
    for (int ish = 0; ish < nshell[iat]; ish++) {
      for (int jsh = 0; jsh < ish; jsh++) {
        gam = hubbard[ish][jsh][izp][izp];
        amat[ii + jsh][ii + ish] += gam;
        amat[ii + ish][ii + jsh] += gam;
      }
      amat[ii + ish][ii + ish] += hubbard[ish][ish][izp][izp];
    }
  }
}

__device__
void es2_update(
  const effective_coulomb &self,
  const structure_type &mol, 
  coulomb_cache &cache)
{
  for (int i = 0; i < MAX_NSH; i++)
  {
    for (int j = 0; j < MAX_NSH; j++)
    {
      cache.amat[i][j] = 0;
    }
  }
  /* if (any(mol%periodic)) then */
      /* periodic not supported */
  /* else */
  /*  call get_amat_0d(mol, self%nshell, self%offset, self%hubbard, self%gexp, amat) */
  get_amat_0d(mol, self.nshell, self.offset, self.hubbard, self.gexp, cache.amat);
}

/*

subroutine ncoord_dexp(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), sigma(3, 3), cutoff2

   cn(:) = 0.0_wp
   dcndr(:, :, :) = 0.0_wp
   dcndL(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   ! $omp parallel do schedule(runtime) default(none) &
   ! $omp reduction(+:cn, dcndr, dcndL) shared(mol, trans, cutoff2, rcov) &
   ! $omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf, countd, sigma)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = exp_count(ka, r1, rc) * exp_count(kb, r1, rc + r_shift)
            countd = (dexp_count(ka, r1, rc) * exp_count(kb, r1, rc + r_shift) &
               & + exp_count(ka, r1, rc) * dexp_count(kb, r1, rc + r_shift)) * rij/r1

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

            dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
            dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
            dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
            dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

            sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

            dcndL(:, :, iat) = dcndL(:, :, iat) + sigma
            if (iat /= jat) then
               dcndL(:, :, jat) = dcndL(:, :, jat) + sigma
            end if

         end do
      end do
   end do

end subroutine ncoord_dexp
*/

__device__
void ncoord_dexp(
  const structure_type &mol,
  /*const float (&trans)[3][1],*/
  const float cutoff,
  const float (&rcov)[MAX_NELEM],
  float (&cn)[MAX_NAT],
  float (&dcndr)[MAX_NAT][MAX_NAT][3],
  float (&dcndL)[MAX_NAT][3][3]
  ) {
    float cutoff2 = cutoff * cutoff;
    float rij[3], countf, countd[3], sigma[3][3];
    int izp, jzp;

    for (int iat = 0; iat < mol.nat; iat++) {
      izp = mol.id[iat];
      cn[iat] = 0.0f;
      for (int jat = 0; jat <= iat; jat++) {
        jzp = mol.id[jat];
        /* for (int itr = 0; itr < MAX_TRANS; itr++) { */
          for (int k = 0; k < 3; k++) {
            rij[k] = mol.xyz[iat][k] - (mol.xyz[jat][k]);
          }
          float r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
          if (r2 > cutoff2 || r2 < 1.0e-12f) continue;

          float r1 = sqrtf(r2);
          float rc = rcov[izp] + rcov[jzp];

          countf = expf(-KA * (r1 - rc)) * expf(-KB * (r1 - (rc + R_SHIFT)));
          for (int k = 0; k < 3; k++) {
            countd[k] = (-KA * expf(-KA * (r1 - rc)) * expf(-KB * (r1 - (rc + R_SHIFT))) +
                         -KB * expf(-KA * (r1 - rc)) * expf(-KB * (r1 - (rc + R_SHIFT)))) *
                        rij[k] / r1;
          }

          cn[iat] += countf;
          if (iat != jat) {
            cn[jat] += countf;
          }

          for (int k = 0; k < 3; k++) {
            dcndr[iat][iat][k] += countd[k];
            dcndr[jat][jat][k] -= countd[k];
            dcndr[iat][jat][k] += countd[k];
            dcndr[jat][iat][k] -= countd[k];
          }

          for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
              sigma[k][l] = countd[k] * rij[l];
            }
          }

          for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
              dcndL[iat][k][l] += sigma[k][l];
              if (iat != jat) {
                dcndL[jat][k][l] += sigma[k][l];
              }
            }
          }
        /* } */
      }
    }
  }

// __device__
__device__
void get_coordination_number(
  const structure_type &mol, 
  /* float (&trans)[3][1], */
  const float cutoff, 
  const float (&rcov)[MAX_NELEM],
  
  float (&cn)[MAX_NAT], 
  float (&dcndr)[MAX_NAT][MAX_NAT][3], 
  float (&dcndL)[MAX_NAT][3][3]
) {
  /* if (present(dcndr) .and. present(dcndL)) then */
    ncoord_dexp(mol, /*trans,*/ cutoff, rcov, cn, dcndr, dcndL);
  /* else */
  /* Unimplemented */
}


__device__
void ncoord_get_cn(
  const gfn_ncoord_type &self,
  const structure_type &mol,
  float (&cn)[MAX_NAT],
  float (&dcndr)[MAX_NAT][MAX_NAT][3],
  float (&dcndL)[MAX_NAT][3][3]
) {
  /* not implemented */
  /* call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, lattr) */

  /* call get_coordination_number(mol, lattr, self%cutoff, self%rcov, cn, dcndr, dcndL) */
  get_coordination_number(mol, /*lattr,*/ self.cutoff, self.rcov, cn, dcndr, dcndL);
}

/* subroutine get_mrad(mol, shift, kexp, rmax, rad, valence_cn, cn, mrad, dmrdcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Shift for the generation of the multipolar damping radii
   real(wp), intent(in) :: shift
   !> Exponent for the generation of the multipolar damping radii
   real(wp), intent(in) :: kexp
   !> Maximum radius for the multipolar damping radii
   real(wp), intent(in) :: rmax
   !> Base radii for the multipolar damping radii
   real(wp), intent(in) :: rad(:)
   !> Valence coordination number
   real(wp), intent(in) :: valence_cn(:)
   !> Coordination numbers for all atoms
   real(wp), intent(in) :: cn(:)
   !> Multipole damping radii for all atoms
   real(wp), intent(out) :: mrad(:)
   !> Derivative of multipole damping radii with repect to the coordination numbers
   real(wp), intent(out) :: dmrdcn(:)

   integer :: iat, izp
   real(wp) :: arg, t1, t2
   ...
*/

__device__
void get_mrad(
  const structure_type &mol,
  const float shift,
  const float kexp,
  const float rmax,
  const float (&rad)[MAX_NELEM],
  const float (&valence_cn)[MAX_NELEM],
  const float (&cn)[MAX_NAT],
  float (&mrad)[MAX_NAT],
  float (&dmrdcn)[MAX_NAT]
)
{
  /* integer :: iat, izp
   real(wp) :: arg, t1, t2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      arg = cn(iat) - valence_cn(izp) - shift
      t1 = exp(-kexp*arg)
      t2 = (rmax - rad(izp)) / (1.0_wp + t1)
      mrad(iat) = rad(izp) + t2
      dmrdcn(iat) = -t2 * kexp * t1 / (1.0_wp + t1)
   end do
   */
  for (int iat = 0; iat < mol.nat; iat++) {
    int izp = mol.id[iat];
    float arg, t1, t2;
    arg = cn[iat] - valence_cn[izp] - shift;
    t1 = expf(-kexp * arg);
    t2 = (rmax - rad[izp]) / (1.0f + t1);
    mrad[iat] = rad[izp] + t2;
    dmrdcn[iat] = -t2 * kexp * t1 / (1.0f + t1);
  }
}

__device__
void aes2_update(
  const damped_multipole &self,
  const structure_type &mol,
  coulomb_cache &ptr
) {
  /* allocs are done within struct itsef */
  /* if (allocated(self%ncoord)) then */
  ncoord_get_cn(self.ncoord, mol, ptr.cn, ptr.dcndr, ptr.dcndL);

  /* get_mrad(mol, self%shift, self%kexp, self%rmax, self%rad, self%valence_cn, &
      & ptr%cn, ptr%mrad, ptr%dmrdcn) */
  get_mrad(mol, self.shift, self.kexp, self.rmax, self.rad, self.valence_cn, 
        ptr.cn, ptr.mrad, ptr.dmrdcn);
}

__device__
void update(
  const tb_coulomb &self, 
  const structure_type &mol,
  coulomb_cache &cache
)
{
  /* if (allocated(self%es2)) then
      call self%es2%update(mol, cache)
   end if */
  es2_update(self.es2, mol, cache);

  /*if (allocated(self%aes2)) then
      call self%aes2%update(mol, cache)
   end if*/
   aes2_update(self.aes2, mol, cache);
}