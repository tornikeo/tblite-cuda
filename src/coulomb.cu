
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
void get_coulomb_matrix(
  const effective_coulomb &self,
  const structure_type &mol, 
  coulomb_cache &cache, 
  float (&amat)[MAX_NSH][MAX_NSH])
{
  for (int i = 0; i < MAX_NSH; i++)
  {
    for (int j = 0; j < MAX_NSH; j++)
    {
      amat[i][j] = 0;
    }
  }
  /* if (any(mol%periodic)) then */
      /* periodic not supported */
  /* else */
  /*  call get_amat_0d(mol, self%nshell, self%offset, self%hubbard, self%gexp, amat) */
  get_amat_0d(mol, self.nshell, self.offset, self.hubbard, self.gexp, amat);
}

__device__
void update(
  const tb_coulomb &self, 
  const structure_type &mol,
  coulomb_cache &cache
)
{
  
}