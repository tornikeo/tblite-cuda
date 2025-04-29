
#include "coulomb.h"
#include "blas/level2.h"
#include <cassert>

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
void effective_coulomb::update(
  // const  &self,
  const structure_type &mol, 
  coulomb_cache &cache) const
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
  get_amat_0d(mol, nshell, offset, hubbard, gexp, cache.amat);
}


__device__ void effective_coulomb::get_potential(
  const structure_type &mol, 
  coulomb_cache &cache, 
  const wavefunction_type &wfn, 
  potential_type &pot) const
{
  symv(MAX_NSH, &cache.amat[0][0], wfn.qsh[1], pot.vsh[1], /*beta=*/1.0);
}

/*!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, ii, ish
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   call symv(ptr%amat, wfn%qsh(:, 1), ptr%vvec, alpha=0.5_wp)
   do iat = 1, mol%nat
      ii = self%offset(iat)
      do ish = 1, self%nshell(iat)
         energies(iat) = energies(iat) + ptr%vvec(ii+ish) * wfn%qsh(ii+ish, 1)
      end do
   end do
end subroutine get_energy
*/
__device__ void effective_coulomb::get_energy(const structure_type &mol, coulomb_cache &cache, const wavefunction_type &wfn, float (&energies)[MAX_NAT]) const
{

  float vvec[MAX_NSH] = {0.0f};
  symv(MAX_NSH, &cache.amat[0][0], wfn.qsh[0], vvec, /*beta=*/0.5f);

  for (int iat = 0; iat < mol.nat; iat++) {
    int ii = offset[iat];
    for (int ish = 0; ish < nshell[iat]; ish++) {
      energies[iat] += vvec[ii + ish] * wfn.qsh[0][ii + ish];
    }
  }
}

__device__
void damped_multipole::update(
  // const damped_multipole &self,
  const structure_type &mol,
  coulomb_cache &ptr
) const {
  /* allocs are done within struct itsef */
  /* if (allocated(self%ncoord)) then */
  ncoord_get_cn(ncoord, mol, ptr.cn, ptr.dcndr, ptr.dcndL);

  /* get_mrad(mol, self%shift, self%kexp, self%rmax, self%rad, self%valence_cn, &
      & ptr%cn, ptr%mrad, ptr%dmrdcn) */
  get_mrad(mol, shift, kexp, rmax, rad, valence_cn, 
        ptr.cn, ptr.mrad, ptr.dmrdcn);

  /* call get_multipole_matrix(self, mol, ptr, ptr%amat_sd, ptr%amat_dd, ptr%amat_sq) */
  get_multipole_matrix(mol, ptr, ptr.amat_sd, ptr.amat_dd, ptr.amat_sq);
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
/* subroutine get_multipole_matrix_0d(mol, rad, kdmp3, kdmp5, amat_sd, amat_dd, amat_sq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole damping radii for all atoms
   real(wp), intent(in) :: rad(:)
   !> Damping function for inverse quadratic contributions
   real(wp), intent(in) :: kdmp3
   !> Damping function for inverse cubic contributions
   real(wp), intent(in) :: kdmp5
   !> Interation matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interation matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interation matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)
  ...
end subroutine get_multipole_matrix_0d
*/

__device__
void get_multipole_matrix_0d(
  const structure_type &mol,
  const float (&rad)[MAX_NAT],
  const float kdmp3,
  const float kdmp5,
  float (&amat_sd)[MAX_NAT][MAX_NAT][3],
  float (&amat_dd)[MAX_NAT][3][MAX_NAT][3],
  float (&amat_sq)[MAX_NAT][MAX_NAT][6]
)
{
  /*
   integer :: iat, jat
   real(wp) :: r1, vec(3), g1, g3, g5, fdmp3, fdmp5, tc(6), rr

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         g1 = 1.0_wp / r1
         g3 = g1 * g1 * g1
         g5 = g3 * g1 * g1

         rr = 0.5_wp * (rad(jat) + rad(iat)) * g1
         fdmp3 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp3)
         fdmp5 = 1.0_wp / (1.0_wp + 6.0_wp * rr**kdmp5)

         amat_sd(:, jat, iat) = amat_sd(:, jat, iat) + vec * g3 * fdmp3
         amat_dd(:, jat, :, iat) = amat_dd(:, jat, :, iat) &
            & + unity * g3*fdmp5 - spread(vec, 1, 3) * spread(vec, 2, 3) * 3*g5*fdmp5
         tc(2) = 2*vec(1)*vec(2)*g5*fdmp5
         tc(4) = 2*vec(1)*vec(3)*g5*fdmp5
         tc(5) = 2*vec(2)*vec(3)*g5*fdmp5
         tc(1) = vec(1)*vec(1)*g5*fdmp5
         tc(3) = vec(2)*vec(2)*g5*fdmp5
         tc(6) = vec(3)*vec(3)*g5*fdmp5
         amat_sq(:, jat, iat) = amat_sq(:, jat, iat) + tc
      end do
   end do */

  for (int iat = 0; iat < mol.nat; iat++) {
    for (int jat = 0; jat < mol.nat; jat++) {
      if (iat == jat) continue;

      float vec[3], r1, g1, g3, g5, rr, fdmp3, fdmp5, tc[6];
      for (int k = 0; k < 3; k++) {
        vec[k] = mol.xyz[iat][k] - mol.xyz[jat][k];
      }

      r1 = sqrtf(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      g1 = 1.0f / r1;
      g3 = g1 * g1 * g1;
      g5 = g3 * g1 * g1;

      rr = 0.5f * (rad[jat] + rad[iat]) * g1;
      fdmp3 = 1.0f / (1.0f + 6.0f * powf(rr, kdmp3));
      fdmp5 = 1.0f / (1.0f + 6.0f * powf(rr, kdmp5));

      for (int k = 0; k < 3; k++) {
        amat_sd[iat][jat][k] += vec[k] * g3 * fdmp3;
      }

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          amat_dd[iat][k][jat][l] += (k == l ? g3 * fdmp5 : 0.0f) -
                                     3.0f * vec[k] * vec[l] * g5 * fdmp5;
        }
      }

      tc[0] = vec[0] * vec[0] * g5 * fdmp5;
      tc[1] = 2.0f * vec[0] * vec[1] * g5 * fdmp5;
      tc[2] = vec[1] * vec[1] * g5 * fdmp5;
      tc[3] = 2.0f * vec[0] * vec[2] * g5 * fdmp5;
      tc[4] = 2.0f * vec[1] * vec[2] * g5 * fdmp5;
      tc[5] = vec[2] * vec[2] * g5 * fdmp5;

      for (int k = 0; k < 6; k++) {
        amat_sq[iat][jat][k] += tc[k];
      }
    }
  }
}


/* subroutine get_multipole_matrix(self, mol, cache, amat_sd, amat_dd, amat_sq)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Interation matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interation matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interation matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)
...
*/
__device__
void damped_multipole::get_multipole_matrix(
  // const damped_multipole &self,
  const structure_type &mol,
  coulomb_cache &cache,
  float (&amat_sd)[MAX_NAT][MAX_NAT][3],
  float (&amat_dd)[MAX_NAT][3][MAX_NAT][3],
  float (&amat_sq)[MAX_NAT][MAX_NAT][6]
) const {
  // amat_sd

  // amat_sd(:, :, :) = 0.0_wp
  // amat_dd(:, :, :, :) = 0.0_wp
  // amat_sq(:, :, :) = 0.0_wp
  for (int i = 0; i < MAX_NAT; i++) {
    for (int j = 0; j < MAX_NAT; j++) {
      for (int k = 0; k < 3; k++) {
        amat_sd[i][j][k] = 0.0f;
        amat_sq[i][j][k] = 0.0f;
        for (int l = 0; l < 3; l++) {
          amat_dd[i][j][k][l] = 0.0f;
        }
      }
    }
  }

  /* if (any(mol%periodic)) then */
  /* Not implemented */
  /* call get_multipole_matrix_3d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & cache%wsc, cache%alpha, amat_sd, amat_dd, amat_sq) */
  /* else */
  /* call get_multipole_matrix_0d(mol, cache%mrad, self%kdmp3, self%kdmp5, &
         & amat_sd, amat_dd, amat_sq) */
  get_multipole_matrix_0d(mol, cache.mrad, kdmp3, kdmp5, amat_sd, amat_dd, amat_sq);
}

/*
!> Get multipolar anisotropic potential contribution
subroutine get_kernel_potential(mol, kernel, mpat, vm)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Potential shoft on atomic multipole moment
   real(wp), intent(inout) :: vm(:, :)

   integer :: iat, izp
   real(wp) :: mpscale(size(mpat, 1)  mpat, vm
   )

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      vm(:, iat) = vm(:, iat) + 2*kernel(izp) * mpat(:, iat) * mpscale
   end do
end subroutine get_kernel_potential
  
*/

template <int D>
__device__
void get_kernel_potential(
  const structure_type &mol,
  const float (&kernel)[MAX_NELEM], 
  const float (&mpat)[MAX_NAT][D], 
  float (&vm)[MAX_NAT][D]
)
{
  float mpscale[D];
  for (int i = 0; i < D; i++) {
    mpscale[i] = 1.0f;
  }
  if (D == 6) {
    mpscale[1] = 2.0f;
    mpscale[3] = 2.0f;
    mpscale[4] = 2.0f;
  }

  for (int iat = 0; iat < mol.nat; iat++) {
    int izp = mol.id[iat];
    for (int d = 0; d < D; d++) {
      vm[iat][d] += 2.0f * kernel[izp] * mpat[iat][d] * mpscale[d];
    }
  }
}


/* 
__device__ void gemv312(const float* A, const float* x, float* y, 
  size_t dim1, size_t dim2, size_t dim3, 
  float alpha = 1.0f, float beta = 0.0f, 
  bool transpose = false)
  {...}
*/
/*
__device__ 
void gemv(const float* A, const float* x, float* y, 
  size_t rows, size_t cols, 
  float alpha = 1.0f, 
  float beta = 0.0f, 
  bool transpose = false) {...}
  */

/*
!> Get multipolar anisotropic exchange-correlation kernel
subroutine get_kernel_energy(mol, kernel, mpat, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp
   real(wp) :: mpt(size(mpat, 1)), mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      mpt(:) = mpat(:, iat) * mpscale
      energies(iat) = energies(iat) + kernel(izp) * dot_product(mpt, mpat(:, iat))
   end do
end subroutine get_kernel_energy
*/


__device__ 
void damped_multipole::get_potential(const structure_type &mol, coulomb_cache &cache, const wavefunction_type &wfn, potential_type &pot) const
{
  /* 
   call gemv(ptr%amat_sd, wfn%qat(:, 1), pot%vdp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sd, wfn%dpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)

   call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call get_kernel_potential(mol, self%dkernel, wfn%dpat(:, :, 1), pot%vdp(:, :, 1))
   call get_kernel_potential(mol, self%qkernel, wfn%qpat(:, :, 1), pot%vqp(:, :, 1))
  */
//  gemv312(&cache.amat_sd[0][0][0], &wfn.qat[0], &pot.vdp[0][0], )
  /* TODO: Go from here */
  /* TODO: Also, you need to do testing - remember the copy-over-plans */
  /* You can copy over data from fortran via printing arrs in fortran */
  gemv(cache.amat_sd, wfn.qat[0], pot.vdp[0],  1.0f, 1.0f, false);
  gemv(cache.amat_sd, wfn.dpat[0], pot.vat[0], 1.0f, 1.0f, true);

  /*     vdp := amat_dd @ dpat    */
  /* 3 9 3 9 x 3 9 -> 3 9 */
  gemv(cache.amat_dd, wfn.dpat[0][0][0], &pot.vdp[0][0][0], mol.nat, 3, mol.nat, 1.0f, 1.0f, false);

  // call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
  // call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")
  gemv(cache.amat_sq, wfn.qat, pot.vqp[0][0][0], mol.nat, mol.nat * 6, 1.0f, 1.0f, false);
  gemv(cache.amat_sq, wfn.qpat, pot.vat[0][0], mol.nat, mol.nat * 6, 1.0f, 1.0f, true);
  
  get_kernel_potential(mol, dkernel, wfn.dpat[0], pot.vdp[0]);
  get_kernel_potential(mol, qkernel, wfn.qpat[0], pot.vqp[0]);
}

/*

!> Get anisotropic electrostatic energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat))

   call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)
   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp, alpha=0.5_wp)
   call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)

   energies(:) = energies + sum(wfn%dpat(:, :, 1) * vd, 1) + sum(wfn%qpat(:, :, 1) * vq, 1)

   call get_kernel_energy(mol, self%dkernel, wfn%dpat(:, :, 1), energies)
   call get_kernel_energy(mol, self%qkernel, wfn%qpat(:, :, 1), energies)
end subroutine get_energy*/

template <int D>
__device__
void get_kernel_energy(
  const structure_type &mol,
  const float (&kernel)[MAX_NELEM], 
  const float (&mpat)[MAX_NAT][D], 
  float (&energies)[MAX_NAT]
)
{
  float mpscale[D];
  for (int i = 0; i < D; i++) {
    mpscale[i] = 1.0f;
  }
  if (D == 6) {
    mpscale[1] = 2.0f;
    mpscale[3] = 2.0f;
    mpscale[4] = 2.0f;
  }

  for (int iat = 0; iat < mol.nat; iat++) {
    int izp = mol.id[iat];
    float dot_product = 0.0f;
    for (int d = 0; d < D; d++) {
      dot_product += mpat[iat][d] * mpat[iat][d] * mpscale[d];
    }
    energies[iat] += kernel[izp] * dot_product;
  }
}

__device__ void damped_multipole::get_energy(const structure_type &mol, coulomb_cache &cache, const wavefunction_type &wfn, float (&energies)[MAX_NAT]) const
{
  float vs[MAX_NAT] {0};
  float vd[MAX_NAT][3] {0};
  float vq[MAX_NAT][6] {0};

  // gemv312(&cache.amat_sd[0][0][0], &wfn.qat[0][0], &pot.vdp[0][0][0], mol.nat, mol.nat * 3, 1.0f, 1.0f, false);
  // gemv321(&cache.amat_sd[0][0][0], &wfn.dpat[0][0][0], &pot.vat[0][0], mol.nat, mol.nat * 3, 1.0f, 1.0f, true);

  // /*     vdp := amat_dd @ dpat    */
  // /* 3 9 3 9 x 3 9 -> 3 9 */
  // gemv422(&cache.amat_dd[0][0][0][0], &wfn.dpat[0][0][0], &pot.vdp[0][0][0], mol.nat, 3, mol.nat, 1.0f, 1.0f, false);

  // // call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
  // // call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")
  // gemv312(&cache.amat_sq[0][0][0], &wfn.qat[0][0], &pot.vqp[0][0][0], mol.nat, mol.nat * 6, 1.0f, 1.0f, false);
  // gemv321(&cache.amat_sq[0][0][0], &wfn.qpat[0][0][0], &pot.vat[0][0], mol.nat, mol.nat * 6, 1.0f, 1.0f, true);
  // 9x9x3 @ 9 + 9x3
  gemv(&cache.amat_sd[0][0][0], &wfn.qat[0][0], &vd[0][0], mol.nat, mol.nat * 3, 1.0f, 0.0f, false);
  // gemv422(&cache.amat_dd[0][0][0][0], &wfn.dpat[0][0][0], &vd[0][0], mol.nat, 3, mol.nat, 0.5f, 1.0f, false);
  // gemv312(&cache.amat_sq[0][0][0], &wfn.qat[0][0], &vq[0][0], mol.nat, mol.nat * 6, 1.0f, 0.0f, false);

  // for (int iat = 0; iat < mol.nat; iat++) {
  //   for (int d = 0; d < 3; d++) {
  //     energies[iat] += wfn.dpat[iat][d][0] * vd[iat][d];
  //   }
  //   for (int q = 0; q < 6; q++) {
  //     energies[iat] += wfn.qpat[iat][q][0] * vq[iat][q];
  //   }
  // }

  // get_kernel_energy(mol, dkernel, wfn.dpat[0], energies);
  // get_kernel_energy(mol, qkernel, wfn.qpat[0], energies);
}

__device__ 
void tb_coulomb::get_potential(
  const structure_type &mol, 
  coulomb_cache &cache,
  const wavefunction_type &wfn,
  potential_type &pot
) const
{
  es2.get_potential(mol, cache, wfn, pot);

  aes2.get_potential(mol, cache, wfn, pot);

  es3.get_potential(mol, cache, wfn, pot);
  /*if (allocated(self%es3)) then
      call self%es3%get_potential(mol, cache, wfn, pot) // UNUSED
   end if*/
}

__device__ 
void tb_coulomb::get_energy(const structure_type &mol, coulomb_cache &cache, const wavefunction_type &wfn, float (&energies)[MAX_NAT]) const
{
  es2.get_energy(mol, cache, wfn, energies);

  aes2.get_energy(mol, cache, wfn, energies);

  // es3.get_energy(mol, cache, wfn, energies);
}

__device__
void tb_coulomb::update(
  // const tb_coulomb &self, 
  const structure_type &mol,
  coulomb_cache &cache
) const
{
  /* if (allocated(self%es2)) then
      call self%es2%update(mol, cache)
   end if */
  es2.update(mol, cache);

  /*if (allocated(self%aes2)) then
      call self%aes2%update(mol, cache)
   end if*/
   aes2.update(mol, cache);
   /* if (allocated(self%es3)) then
      call self%es3%update(mol, cache)
      Not allocated!
   end if
   */

}


void onsite_thirdorder::get_potential(
  const structure_type &mol, const coulomb_cache &cache, 
  const wavefunction_type &wfn, potential_type &pot) const
{
  if (shell_resolved) {
    for (int iat = 0; iat < mol.nat; iat++) {
      int izp = mol.id[iat];
      int ii = ish_at[iat];
      for (int ish = 0; ish < nsh_at[iat]; ish++) {
        pot.vsh[0][ii + ish] += wfn.qsh[0][ii + ish] * wfn.qsh[0][ii + ish] * hubbard_derivs[izp][ish];
      }
    }
  } else {
    for (int iat = 0; iat < mol.nat; iat++) {
      int izp = mol.id[iat];
      pot.vat[0][iat] += wfn.qat[0][iat] * wfn.qat[0][iat] * hubbard_derivs[izp][0];
    }
  }
}

/*
!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electrostatic container
   class(onsite_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp, ii, ish

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            energies(iat) = energies(iat) &
               & + wfn%qsh(ii+ish, 1)**3 * self%hubbard_derivs(ish, izp) / 3.0_wp
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         energies(iat) = energies(iat) &
            & + wfn%qat(iat, 1)**3 * self%hubbard_derivs(1, izp) / 3.0_wp
      end do
   end if
end subroutine get_energy
*/
__device__ void onsite_thirdorder::get_energy(const structure_type &mol, coulomb_cache &cache, const wavefunction_type &wfn, float (&energies)[MAX_NAT]) const
{
  /*if (shell_resolved) {*/ assert(shell_resolved&&"We only support shell resolved in this function");
  for (int iat = 0; iat < mol.nat; iat++) {
    int izp = mol.id[iat];
    int ii = ish_at[iat];
    for (int ish = 0; ish < nsh_at[iat]; ish++) {
      energies[iat] += powf(wfn.qsh[0][ii + ish], 3) * hubbard_derivs[izp][ish] / 3.0f;
    }
  }
  /*} else {
    for (int iat = 0; iat < mol.nat; iat++) {
      int izp = mol.id[iat];
      energies[iat] += powf(wfn.qat[0][iat], 3) * hubbard_derivs[izp][0] / 3.0f;
    }
  }*/
}
